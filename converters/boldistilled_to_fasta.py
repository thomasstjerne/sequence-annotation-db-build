"""
Convert BOLDistilled FASTA + taxonomy TSV to the normalised FASTA header format.

BOLDistilled FASTA header format:  >sampleID|BIN
Taxonomy TSV is indexed by BIN with columns:
  bin, kingdom, phylum, class, order, family, subfamily, tribe,
  genus, species, subspecies, concordant_rank, discordant_rank

Output header format (pipe-separated, same field order as all other conversion scripts):
  ID | accessionNumber | scientificName | decimalLatitude | decimalLongitude |
  typeStatus | catalogueNumber | identifiedBy | taxonRank | country | locality |
  basisOfRecord | higherClassification | dataset | targetGene

BOLDistilled releases carry the release month in the filename
(BOLDistilled_COI_<Mon><YYYY>_SEQUENCES.fasta), and the endpoint serves whichever
release is current. Point this script at the extracted source directory and it
resolves the newest release itself, so no version needs pinning in datasets.yaml.

Usage:
    python3 converters/boldistilled_to_fasta.py <source_dir> <dataset_shortname> --target-gene coi
    python3 converters/boldistilled_to_fasta.py <fasta_file> <taxonomy_tsv> <dataset_shortname> --target-gene coi

Example:
    python3 converters/boldistilled_to_fasta.py \\
        source-data/boldistilled/source boldistilled_coi --target-gene coi
"""

import argparse
import csv
import re
import unicodedata
from pathlib import Path

OUTPUT_DIR = Path("output/fasta")

TAXON_FIELDS = ["kingdom", "phylum", "class", "order", "family", "genus", "species"]

# BOLDistilled_COI_Apr2026_SEQUENCES.fasta -> version "Apr2026"
RELEASE_RE = re.compile(
    r"^BOLDistilled_(?P<gene>[A-Za-z0-9]+)_(?P<version>(?P<mon>[A-Za-z]{3})(?P<year>\d{4}))_SEQUENCES\.fasta$",
    re.IGNORECASE,
)

MONTHS = {m: i for i, m in enumerate(
    ["jan", "feb", "mar", "apr", "may", "jun",
     "jul", "aug", "sep", "oct", "nov", "dec"], start=1)}



# '|' separates header fields and '>' starts a FASTA record, so neither may appear
# inside a field value. Source data does contain them: NBDL's identifiedBy lists
# multiple collectors as "Pogonoski | Russell", which silently shifted every later
# field for 28 records until a header-width check caught it.
def strip_delimiters(value: str) -> str:
    return str(value).replace("|", "/").replace(">", "") if value else ""

def release_sort_key(match: re.Match, path: Path) -> tuple:
    """Order releases newest-last. Unparseable months fall back to file mtime."""
    year = int(match.group("year"))
    month = MONTHS.get(match.group("mon").lower())
    if month is None:
        return (0, 0, path.stat().st_mtime)
    return (year, month, path.stat().st_mtime)


def discover_release(directory: Path) -> tuple[Path, Path, str]:
    """Find the newest BOLDistilled release in *directory*.

    Returns (sequences_fasta, taxonomy_tsv, version). Raises if none is usable.
    """
    candidates = []
    for fasta in directory.iterdir():
        m = RELEASE_RE.match(fasta.name)
        if not m:
            continue
        tsv = fasta.with_name(fasta.name.replace("_SEQUENCES.fasta", "_TAXONOMY.tsv"))
        if not tsv.exists():
            print(f"  ! {fasta.name}: no matching _TAXONOMY.tsv, skipping")
            continue
        candidates.append((release_sort_key(m, fasta), fasta, tsv, m.group("version")))

    if not candidates:
        raise FileNotFoundError(
            f"No BOLDistilled_*_SEQUENCES.fasta with a matching _TAXONOMY.tsv found in {directory}"
        )

    candidates.sort(key=lambda c: c[0])
    if len(candidates) > 1:
        others = ", ".join(c[3] for c in candidates[:-1])
        print(f"  {len(candidates)} releases present ({others}, {candidates[-1][3]}) — using newest")

    _, fasta, tsv, version = candidates[-1]
    return fasta, tsv, version


def sanitize(value) -> str:
    if not value or value in ("", "None", "nan"):
        return ""
    normalized = unicodedata.normalize("NFD", str(value))
    ascii_only = normalized.encode("ascii", "ignore").decode("ascii")
    ascii_only = strip_delimiters(ascii_only)
    return re.sub(r"\s+", "_", ascii_only).strip("_")


def get_scientific_name(row: dict, taxon_rank: str) -> str:
    """Return the name at the concordant rank, falling back down the hierarchy."""
    if taxon_rank and taxon_rank != "None":
        val = row.get(taxon_rank, "")
        if val and val != "None":
            return sanitize(val)
    # Fall back to lowest non-empty standard rank
    for field in reversed(TAXON_FIELDS):
        val = row.get(field, "")
        if val and val != "None":
            return sanitize(val)
    return ""


def get_higher_classification(row: dict) -> str:
    parts = []
    for field in TAXON_FIELDS:
        val = row.get(field, "")
        if val and val != "None":
            parts.append(sanitize(val))
    return ";".join(parts)


def load_taxonomy(tsv_path: Path) -> dict[str, dict]:
    taxonomy: dict[str, dict] = {}
    with open(tsv_path, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            bin_id = row.get("bin", "").strip()
            if bin_id:
                taxonomy[bin_id] = row
    return taxonomy


def build_header(sample_id: str, bin_id: str, tax: dict, dataset: str, target_gene: str) -> str:
    taxon_rank    = tax.get("concordant_rank", "")
    if taxon_rank == "None":
        taxon_rank = ""
    scientific_name        = get_scientific_name(tax, taxon_rank)
    higher_classification  = get_higher_classification(tax)

    def rank(field: str) -> str:
        val = tax.get(field, "")
        return sanitize(val) if val and val != "None" else ""

    fields = [
        sample_id,
        sample_id,
        scientific_name,
        "",                          # decimalLatitude
        "",                          # decimalLongitude
        "",                          # typeStatus
        "",                          # catalogueNumber
        "",                          # identifiedBy
        taxon_rank,
        "",                          # country
        "",                          # locality
        "",                          # basisOfRecord
        higher_classification,
        dataset,
        target_gene,
        "",                          # domain
        rank("kingdom"),
        rank("phylum"),
        rank("class"),
        rank("order"),
        rank("family"),
        rank("genus"),
        rank("species"),
    ]
    return "|".join(strip_delimiters(f) for f in fields)


def main():
    parser = argparse.ArgumentParser(description="Convert BOLDistilled to normalised FASTA format")
    parser.add_argument(
        "inputs", nargs="+",
        help="Either <source_dir> <dataset>, or <fasta_file> <taxonomy_tsv> <dataset>",
    )
    parser.add_argument("--target-gene", required=True, help="Target gene label (e.g. coi)")
    parser.add_argument("--output-dir",  default=None, help="Directory to write output FASTA (default: output/fasta)")
    args = parser.parse_args()

    version = None
    if len(args.inputs) == 2:
        source_dir, dataset = Path(args.inputs[0]), args.inputs[1]
        if not source_dir.is_dir():
            parser.error(f"Not a directory: {source_dir} (pass a source dir, or fasta + tsv + dataset)")
        print(f"Resolving BOLDistilled release in {source_dir} …")
        fasta_path, tsv_path, version = discover_release(source_dir)
        print(f"  release {version}: {fasta_path.name}")
    elif len(args.inputs) == 3:
        fasta_path, tsv_path, dataset = Path(args.inputs[0]), Path(args.inputs[1]), args.inputs[2]
        m = RELEASE_RE.match(fasta_path.name)
        version = m.group("version") if m else None
    else:
        parser.error("expected 2 arguments (source_dir dataset) or 3 (fasta tsv dataset)")

    OUTPUT_DIR  = Path(args.output_dir) if args.output_dir else Path("output/fasta")
    target_gene = args.target_gene.lower()
    output_path = OUTPUT_DIR / f"{dataset}.fasta"

    for p in (fasta_path, tsv_path):
        if not p.exists():
            raise FileNotFoundError(f"File not found: {p}")

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    print("Loading taxonomy …")
    taxonomy = load_taxonomy(tsv_path)
    print(f"  {len(taxonomy):,} BIN records loaded")

    print(f"Streaming {fasta_path.name} …")
    written = missing = 0

    with open(fasta_path, encoding="utf-8") as inp, \
         open(output_path, "w", encoding="utf-8") as out:

        sample_id = bin_id = ""
        seq_lines = []

        def flush():
            nonlocal written, missing
            if not sample_id:
                return
            tax = taxonomy.get(bin_id)
            if not tax:
                nonlocal missing
                missing += 1
                return
            header = build_header(sample_id, bin_id, tax, dataset, target_gene)
            out.write(f">{header}\n{''.join(seq_lines).upper()}\n")
            written += 1

        for line in inp:
            line = line.rstrip("\n")
            if line.startswith(">"):
                flush()
                seq_lines = []
                raw = line[1:].split()[0]          # e.g. BLPAA17045-20|BOLD:AAA0017
                parts = raw.split("|", 1)
                sample_id = parts[0]
                bin_id    = parts[1] if len(parts) > 1 else ""
            else:
                seq_lines.append(line.strip())

        flush()

    print(f"\nDone — {written:,} sequences written to {output_path}")
    if version:
        print(f"  source release: BOLDistilled {version}")
        Path(f"{output_path}.version").write_text(f"boldistilled\t{version}\n")
    if missing:
        print(f"  {missing:,} skipped (BIN not found in taxonomy TSV)")


if __name__ == "__main__":
    main()
