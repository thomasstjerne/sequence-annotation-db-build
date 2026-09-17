"""
Download matK (or another gene) reference sequences from NCBI Entrez and write
them in the normalised 23-field header format.

Unlike the file-based converters, this script fetches directly from NCBI:
  1. esearch (with history) for the gene query  → WebEnv / QueryKey / count
  2. efetch FASTA in batches                    → sequences + accessions
  3. esummary in batches                        → accession → taxid
  4. efetch db=taxonomy for the unique taxids   → ranked lineage per taxid
  5. join + write the normalised FASTA

Raw downloads are cached under --raw-dir so re-runs skip the network unless
--force is given.

API key (fail-tolerant): if NCBI_API_KEY is set (env or .env at repo root) it is
used and the rate limit is 10 req/s; otherwise the script still runs at 3 req/s
and prints a notice recommending you configure one. NCBI_EMAIL / NCBI_TOOL are
sent when present (NCBI etiquette for bulk use).

Output header format (pipe-separated, same field order as dwc_to_fasta.py):
  ID | accessionNumber | scientificName | decimalLatitude | decimalLongitude |
  typeStatus | catalogueNumber | identifiedBy | taxonRank | country | locality |
  basisOfRecord | higherClassification | dataset | targetGene |
  domain | kingdom | phylum | class | order | family | genus | species

Usage:
    python3 converters/ncbi_matk_to_fasta.py ncbi_matk --target-gene matK \\
        --raw-dir source-data/ncbi_matk

    # quick test against a small slice:
    python3 converters/ncbi_matk_to_fasta.py ncbi_matk --target-gene matK \\
        --raw-dir /tmp/ncbi_matk_test --limit 200
"""

import argparse
import sys
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

DEFAULT_QUERY = "matK[Gene Name] AND Viridiplantae[Organism] AND 250:2000[SLEN]"

# NCBI taxonomy rank → our header field.
# NCBI now ranks Eukaryota/Bacteria/Archaea as "domain" (formerly "superkingdom");
# accept both so the domain field is populated.
RANK_TO_FIELD = {
    "domain": "domain",
    "superkingdom": "domain",
    "kingdom": "kingdom",
    "phylum": "phylum",
    "class": "class",
    "order": "order",
    "family": "family",
    "genus": "genus",
    "species": "species",
}


# ── .env loading (no dependency) ──────────────────────────────────────────────


# '|' separates header fields and '>' starts a FASTA record, so neither may appear
# inside a field value. Source data does contain them: NBDL's identifiedBy lists
# multiple collectors as "Pogonoski | Russell", which silently shifted every later
# field for 28 records until a header-width check caught it.
def strip_delimiters(value: str) -> str:
    return str(value).replace("|", "/").replace(">", "") if value else ""

def load_dotenv():
    """Populate os.environ from a .env at the repo root, without overriding real env."""
    import os
    repo_root = Path(__file__).resolve().parent.parent
    env_path = repo_root / ".env"
    if not env_path.exists():
        return
    for line in env_path.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        key, val = line.split("=", 1)
        os.environ.setdefault(key.strip(), val.strip())


# ── Entrez helpers ────────────────────────────────────────────────────────────

class Entrez:
    def __init__(self):
        import os
        self.api_key = os.getenv("NCBI_API_KEY") or ""
        self.email = os.getenv("NCBI_EMAIL") or ""
        self.tool = os.getenv("NCBI_TOOL") or "sequence-annotation-db-build"
        # NCBI allows 10 req/s with a key, 3 req/s without.
        self.delay = 0.11 if self.api_key else 0.34
        self._last = 0.0
        if self.api_key:
            print(f"NCBI API key detected — requesting at ~10 req/s (tool={self.tool}).")
        else:
            print(
                "NOTICE: no NCBI_API_KEY configured — running at the slower anonymous "
                "rate (~3 req/s).\n"
                "        It is recommended to configure one: copy .env.example to .env "
                "and set NCBI_API_KEY\n"
                "        (free, from https://www.ncbi.nlm.nih.gov/account/ → API Key "
                "Management). NCBI_EMAIL is also recommended.",
                file=sys.stderr,
            )

    def _common(self):
        p = {"tool": self.tool}
        if self.api_key:
            p["api_key"] = self.api_key
        if self.email:
            p["email"] = self.email
        return p

    def _throttle(self):
        wait = self.delay - (time.monotonic() - self._last)
        if wait > 0:
            time.sleep(wait)
        self._last = time.monotonic()

    def get(self, endpoint, params, retries=4):
        params = {**self._common(), **params}
        url = f"{EUTILS}/{endpoint}?{urllib.parse.urlencode(params)}"
        for attempt in range(retries):
            self._throttle()
            try:
                with urllib.request.urlopen(url, timeout=120) as resp:
                    return resp.read().decode("utf-8", "replace")
            except Exception as e:  # noqa: BLE001 — transient HTTP/network errors
                if attempt == retries - 1:
                    raise
                time.sleep(1.5 * (attempt + 1))

    def post(self, endpoint, params, retries=4):
        params = {**self._common(), **params}
        data = urllib.parse.urlencode(params).encode("utf-8")
        url = f"{EUTILS}/{endpoint}"
        for attempt in range(retries):
            self._throttle()
            try:
                with urllib.request.urlopen(url, data=data, timeout=120) as resp:
                    return resp.read().decode("utf-8", "replace")
            except Exception:  # noqa: BLE001
                if attempt == retries - 1:
                    raise
                time.sleep(1.5 * (attempt + 1))


# ── pipeline steps ────────────────────────────────────────────────────────────

def esearch(ez, query):
    xml = ez.get("esearch.fcgi", {"db": "nuccore", "term": query, "usehistory": "y", "retmax": 0})
    root = ET.fromstring(xml)
    count = int(root.findtext("Count", "0"))
    return count, root.findtext("WebEnv", ""), root.findtext("QueryKey", "")


def fetch_fasta(ez, webenv, qk, total, raw_path, batch):
    written = 0
    with open(raw_path, "w", encoding="utf-8") as out:
        for start in range(0, total, batch):
            text = ez.get("efetch.fcgi", {
                "db": "nuccore", "WebEnv": webenv, "query_key": qk,
                "rettype": "fasta", "retmode": "text",
                "retstart": start, "retmax": min(batch, total - start),
            })
            out.write(text if text.endswith("\n") else text + "\n")
            written = min(start + batch, total)
            print(f"  efetch fasta {written:,}/{total:,}", end="\r", flush=True)
    print()
    return raw_path


def fetch_taxid_map(ez, webenv, qk, total, batch):
    """accession.version → taxid via esummary over the search history."""
    acc2taxid = {}
    for start in range(0, total, batch):
        xml = ez.get("esummary.fcgi", {
            "db": "nuccore", "WebEnv": webenv, "query_key": qk,
            "version": "2.0", "retstart": start, "retmax": min(batch, total - start),
        })
        root = ET.fromstring(xml)
        for ds in root.iter("DocumentSummary"):
            acc = ds.findtext("AccessionVersion")
            taxid = ds.findtext("TaxId")
            if acc and taxid:
                acc2taxid[acc] = taxid
        print(f"  esummary taxids {min(start + batch, total):,}/{total:,}", end="\r", flush=True)
    print()
    return acc2taxid


def fetch_lineages(ez, taxids, batch=200):
    """taxid → {field: name} ranked lineage via efetch db=taxonomy."""
    lineages = {}
    taxids = list(taxids)
    for start in range(0, len(taxids), batch):
        chunk = taxids[start:start + batch]
        xml = ez.post("efetch.fcgi", {
            "db": "taxonomy", "id": ",".join(chunk), "retmode": "xml",
        })
        root = ET.fromstring(xml)
        for taxon in root.findall("Taxon"):
            tid = taxon.findtext("TaxId")
            fields = {}
            higher = []
            for node in taxon.findall("./LineageEx/Taxon"):
                name = node.findtext("ScientificName", "")
                rank = node.findtext("Rank", "")
                if name:
                    higher.append(name)
                if rank in RANK_TO_FIELD and name:
                    fields[RANK_TO_FIELD[rank]] = name
            # the terminal taxon itself (usually the species)
            term_name = taxon.findtext("ScientificName", "")
            term_rank = taxon.findtext("Rank", "")
            if term_name:
                higher.append(term_name)
                if term_rank in RANK_TO_FIELD:
                    fields[RANK_TO_FIELD[term_rank]] = term_name
            fields["_scientificName"] = term_name
            fields["_taxonRank"] = term_rank if term_rank != "no rank" else ""
            fields["_higher"] = ";".join(higher)
            lineages[tid] = fields
        print(f"  taxonomy {min(start + batch, len(taxids)):,}/{len(taxids):,}", end="\r", flush=True)
    print()
    return lineages


# ── header building ───────────────────────────────────────────────────────────

def build_header(accession, tax, dataset, target_gene):
    def s(v):
        return v.replace(" ", "_") if v else ""

    tax = tax or {}
    fields = [
        accession,                         # ID
        accession,                         # accessionNumber
        s(tax.get("_scientificName", "")),  # scientificName
        "", "", "", "", "",                # lat, lon, typeStatus, catalogue, identifiedBy
        tax.get("_taxonRank", ""),         # taxonRank
        "", "", "",                        # country, locality, basisOfRecord
        s(tax.get("_higher", "")),         # higherClassification
        dataset,                           # dataset
        target_gene,                       # targetGene (kept verbatim, e.g. matK)
        s(tax.get("domain", "")),
        s(tax.get("kingdom", "")),
        s(tax.get("phylum", "")),
        s(tax.get("class", "")),
        s(tax.get("order", "")),
        s(tax.get("family", "")),
        s(tax.get("genus", "")),
        s(tax.get("species", "")),
    ]
    return "|".join(strip_delimiters(f) for f in fields)


def iter_fasta(path):
    """Yield (accession, [seq_lines]) from a raw NCBI FASTA."""
    acc, lines = None, []
    with open(path, encoding="utf-8") as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if acc:
                    yield acc, lines
                acc = line[1:].split()[0] if len(line) > 1 else ""
                lines = []
            else:
                lines.append(line.strip())
    if acc:
        yield acc, lines


# ── main ──────────────────────────────────────────────────────────────────────

def main():
    p = argparse.ArgumentParser(description="Download a gene from NCBI Entrez → normalised FASTA")
    p.add_argument("dataset", help="Short dataset name for headers/output (e.g. ncbi_matk)")
    p.add_argument("--target-gene", required=True, help="targetGene label written to headers (e.g. matK)")
    p.add_argument("--query", default=DEFAULT_QUERY, help="Entrez nuccore query (default: matK / Viridiplantae)")
    p.add_argument("--output-dir", default=None, help="Directory for output FASTA (default: output/fasta)")
    p.add_argument("--raw-dir", default=None, help="Cache dir for raw downloads (default: source-data/<dataset>)")
    p.add_argument("--batch-size", type=int, default=500, help="efetch/esummary batch size")
    p.add_argument("--limit", type=int, default=None, help="Cap number of records (for testing)")
    p.add_argument("--force", action="store_true", help="Re-download even if cached raw files exist")
    args = p.parse_args()

    load_dotenv()
    ez = Entrez()

    output_dir = Path(args.output_dir) if args.output_dir else Path("output/fasta")
    raw_dir = Path(args.raw_dir) if args.raw_dir else Path("source-data") / args.dataset
    output_dir.mkdir(parents=True, exist_ok=True)
    raw_dir.mkdir(parents=True, exist_ok=True)

    raw_fasta = raw_dir / f"{args.dataset}_raw.fasta"
    taxid_tsv = raw_dir / f"{args.dataset}_acc_taxid.tsv"

    # 1. esearch
    print(f"esearch: {args.query}")
    count, webenv, qk = esearch(ez, args.query)
    total = min(count, args.limit) if args.limit else count
    print(f"  {count:,} records found; fetching {total:,}.")
    if total == 0:
        print("Nothing to fetch — exiting.")
        return

    # 2. efetch fasta (+ cache)
    if args.force or not raw_fasta.exists():
        print("Downloading sequences …")
        fetch_fasta(ez, webenv, qk, total, raw_fasta, args.batch_size)
    else:
        print(f"Using cached sequences: {raw_fasta}")

    # 3. accession → taxid (+ cache)
    if args.force or not taxid_tsv.exists():
        print("Resolving taxids …")
        acc2taxid = fetch_taxid_map(ez, webenv, qk, total, args.batch_size)
        with open(taxid_tsv, "w", encoding="utf-8") as f:
            for acc, tid in acc2taxid.items():
                f.write(f"{acc}\t{tid}\n")
    else:
        print(f"Using cached taxids: {taxid_tsv}")
        acc2taxid = dict(l.rstrip("\n").split("\t") for l in open(taxid_tsv) if "\t" in l)

    # 4. taxid → ranked lineage
    print("Fetching ranked lineages …")
    lineages = fetch_lineages(ez, set(acc2taxid.values()))

    # 5. join + write normalised FASTA
    output_path = output_dir / f"{args.dataset}.fasta"
    written = no_tax = 0
    with open(output_path, "w", encoding="utf-8") as out:
        for acc, seq_lines in iter_fasta(raw_fasta):
            seq = "".join(seq_lines).upper()
            if not seq:
                continue
            taxid = acc2taxid.get(acc)
            tax = lineages.get(taxid) if taxid else None
            if not tax:
                no_tax += 1
            header = build_header(acc, tax, args.dataset, args.target_gene)
            out.write(f">{header}\n{seq}\n")
            written += 1

    print(f"Done — {written:,} sequences written to {output_path}")
    if no_tax:
        print(f"  ({no_tax:,} records had no resolvable taxonomy — written with empty rank fields)")


if __name__ == "__main__":
    main()
