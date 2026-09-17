"""
Convert GTDB SSU FASTA files + taxonomy TSVs to the normalised FASTA header format.

GTDB provides separate files for bacteria (bac120) and archaea (ar53).
Both are auto-discovered from the input directory via glob.

GTDB SSU FASTA header format:
  >GB_GCA_003054575.1~DCZV01000001.1_1 [optional bracketed metadata]
  Genome accession = prefix before '~' (or first whitespace-delimited token if no '~')

GTDB taxonomy TSV format (no header row):
  GB_GCA_003054575.1<TAB>d__Bacteria;p__Proteobacteria;c__...;o__...;f__...;g__...;s__...

Rank prefixes: d__=domain, p__=phylum, c__=class, o__=order,
               f__=family, g__=genus, s__=species

Output header format (pipe-separated, same field order as all other conversion scripts):
  ID | accessionNumber | scientificName | decimalLatitude | decimalLongitude |
  typeStatus | catalogueNumber | identifiedBy | taxonRank | country | locality |
  basisOfRecord | higherClassification | dataset | targetGene

Usage:
    python3 converters/gtdb_to_fasta.py <gtdb_dir> <dataset_shortname> --target-gene 16s

Example:
    python3 converters/gtdb_to_fasta.py source-data/gtdb gtdb_16s --target-gene 16s
"""

import argparse
import gzip
import re
import unicodedata
from pathlib import Path

OUTPUT_DIR = Path("output/fasta")

RANK_PREFIXES = {
    "d__": "domain",
    "p__": "phylum",
    "c__": "class",
    "o__": "order",
    "f__": "family",
    "g__": "genus",
    "s__": "species",
}

RANK_ORDER = ["domain", "phylum", "class", "order", "family", "genus", "species"]

# Placeholder patterns used by GTDB for unresolved nodes
UNRESOLVED = re.compile(r"^[a-z]__$")



# '|' separates header fields and '>' starts a FASTA record, so neither may appear
# inside a field value. Source data does contain them: NBDL's identifiedBy lists
# multiple collectors as "Pogonoski | Russell", which silently shifted every later
# field for 28 records until a header-width check caught it.
def strip_delimiters(value: str) -> str:
    return str(value).replace("|", "/").replace(">", "") if value else ""

def sanitize(value: str) -> str:
    if not value:
        return ""
    normalized = unicodedata.normalize("NFD", value)
    ascii_only = normalized.encode("ascii", "ignore").decode("ascii")
    ascii_only = strip_delimiters(ascii_only)
    return re.sub(r"\s+", "_", ascii_only).strip("_")


def parse_taxonomy(tax_str: str) -> dict[str, str]:
    """Parse 'd__X;p__Y;...' into {rank: name}, skipping empty/unresolved nodes."""
    ranks: dict[str, str] = {}
    for node in tax_str.split(";"):
        node = node.strip()
        for prefix, rank in RANK_PREFIXES.items():
            if node.startswith(prefix):
                name = node[len(prefix):]
                if name and not UNRESOLVED.match(node):
                    ranks[rank] = name
                break
    return ranks


def get_scientific_name_and_rank(ranks: dict[str, str]) -> tuple[str, str]:
    for rank in reversed(RANK_ORDER):
        if rank in ranks:
            return sanitize(ranks[rank]), rank
    return "", ""


def get_higher_classification(ranks: dict[str, str]) -> str:
    return ";".join(sanitize(ranks[r]) for r in RANK_ORDER if r in ranks)


def load_taxonomy(*tsv_paths: Path) -> dict[str, dict[str, str]]:
    """Load one or more taxonomy TSVs (no header) into a genome-accession → ranks dict."""
    taxonomy: dict[str, dict[str, str]] = {}
    for path in tsv_paths:
        with open(path, encoding="utf-8") as f:
            for line in f:
                line = line.rstrip("\n")
                if not line:
                    continue
                parts = line.split("\t", 1)
                if len(parts) < 2:
                    continue
                accession, tax_str = parts[0].strip(), parts[1].strip()
                taxonomy[accession] = parse_taxonomy(tax_str)
    return taxonomy


def genome_accession(fasta_header: str) -> str:
    """Extract genome accession from GTDB SSU FASTA header (part before '~' or first space)."""
    token = fasta_header.lstrip(">").split()[0]
    return token.split("~")[0]


def build_header(seq_id: str, accession: str, ranks: dict[str, str],
                 dataset: str, target_gene: str) -> str:
    scientific_name, taxon_rank = get_scientific_name_and_rank(ranks)
    higher_classification       = get_higher_classification(ranks)

    def r(rank: str) -> str:
        return sanitize(ranks.get(rank, ""))

    fields = [
        seq_id,
        accession,
        scientific_name,
        "",   # decimalLatitude
        "",   # decimalLongitude
        "",   # typeStatus
        "",   # catalogueNumber
        "",   # identifiedBy
        taxon_rank,
        "",   # country
        "",   # locality
        "",   # basisOfRecord
        higher_classification,
        dataset,
        target_gene,
        r("domain"),
        "",           # kingdom (not in GTDB)
        r("phylum"),
        r("class"),
        r("order"),
        r("family"),
        r("genus"),
        r("species"),
    ]
    return "|".join(strip_delimiters(f) for f in fields)


def open_fasta(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return open(path, encoding="utf-8")


def convert_fasta(fasta_path: Path, taxonomy: dict[str, dict[str, str]],
                  dataset: str, target_gene: str, out, stats: dict) -> None:
    accession = seq_id = ""
    seq_lines: list[str] = []

    def flush():
        if not seq_id or not seq_lines:
            return
        ranks = taxonomy.get(accession)
        if not ranks:
            stats["missing"] += 1
            return
        header = build_header(seq_id, accession, ranks, dataset, target_gene)
        out.write(f">{header}\n{''.join(seq_lines).upper()}\n")
        stats["written"] += 1

    with open_fasta(fasta_path) as inp:
        for line in inp:
            line = line.rstrip("\n")
            if line.startswith(">"):
                flush()
                seq_lines = []
                accession = genome_accession(line)
                seq_id    = line.lstrip(">").split()[0]
            else:
                seq_lines.append(line.strip())
        flush()


def find_glob(directory: Path, pattern: str) -> list[Path]:
    return sorted(directory.glob(pattern))


def main():
    parser = argparse.ArgumentParser(description="Convert GTDB SSU FASTA to normalised header format")
    parser.add_argument("gtdb_dir",    help="Directory containing GTDB SSU FASTA and taxonomy files")
    parser.add_argument("dataset",     help="Short dataset name for headers (e.g. gtdb_16s)")
    parser.add_argument("--target-gene", required=True, help="Target gene label (e.g. 16s)")
    parser.add_argument("--output-dir",  default=None,  help="Directory to write output FASTA (default: output/fasta)")
    args = parser.parse_args()

    OUTPUT_DIR  = Path(args.output_dir) if args.output_dir else Path("output/fasta")
    gtdb_dir    = Path(args.gtdb_dir)
    target_gene = args.target_gene.lower()
    output_path = OUTPUT_DIR / f"{args.dataset}.fasta"

    # Auto-discover files — handles both versioned (bac120_ssu_reps_r220.fna.gz)
    # and unversioned (bac120_ssu_reps.fna.gz) filenames
    fasta_paths = (find_glob(gtdb_dir, "bac120_ssu_reps*.fna.gz") +
                   find_glob(gtdb_dir, "bac120_ssu_reps*.fna") +
                   find_glob(gtdb_dir, "ar53_ssu_reps*.fna.gz") +
                   find_glob(gtdb_dir, "ar53_ssu_reps*.fna"))

    tax_paths = (find_glob(gtdb_dir, "bac120_taxonomy*.tsv") +
                 find_glob(gtdb_dir, "ar53_taxonomy*.tsv"))

    if not fasta_paths:
        raise FileNotFoundError(f"No bac120/ar53 SSU FASTA files found in {gtdb_dir}")
    if not tax_paths:
        raise FileNotFoundError(f"No bac120/ar53 taxonomy TSV files found in {gtdb_dir}")

    print(f"FASTA files : {[p.name for p in fasta_paths]}")
    print(f"Taxonomy    : {[p.name for p in tax_paths]}")

    print("Loading taxonomy …")
    taxonomy = load_taxonomy(*tax_paths)
    print(f"  {len(taxonomy):,} genome records loaded")

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    stats = {"written": 0, "missing": 0}

    print(f"Writing {output_path} …")
    with open(output_path, "w", encoding="utf-8") as out:
        for fasta_path in fasta_paths:
            print(f"  Processing {fasta_path.name} …")
            convert_fasta(fasta_path, taxonomy, args.dataset, target_gene, out, stats)

    print(f"\nDone — {stats['written']:,} sequences written to {output_path}")
    if stats["missing"]:
        print(f"  {stats['missing']:,} skipped (genome accession not found in taxonomy)")


if __name__ == "__main__":
    main()
