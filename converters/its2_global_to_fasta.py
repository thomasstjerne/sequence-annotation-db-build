"""
Convert the ITS2 Global database FASTA to the normalised header format.

Source FASTA header format:
  >ID;tax=k:Kingdom,p:Phylum,c:Class,o:Order,f:Family,g:Genus,s:Species_name;

Species names use underscores for spaces (e.g. Aesculus_hippocastanum).

Output header format (pipe-separated, same field order as all other conversion scripts):
  ID | accessionNumber | scientificName | decimalLatitude | decimalLongitude |
  typeStatus | catalogueNumber | identifiedBy | taxonRank | country | locality |
  basisOfRecord | higherClassification | dataset | targetGene |
  domain | kingdom | phylum | class | order | family | genus | species

Reference:
  Quaresma et al. (2024) Semi-automated sequence curation for reliable reference
  datasets in ITS2 vascular plant DNA (meta-)barcoding. Scientific Data.
  https://doi.org/10.1038/s41597-024-02962-5

Usage:
    python3 converters/its2_global_to_fasta.py <fasta_file> <dataset_shortname> --target-gene ITS2

Example:
    python3 converters/its2_global_to_fasta.py \\
        source-data/its2_global/its2.global.2023-01-17.curated.tax.mc.add.fa \\
        its2_global --target-gene ITS2
"""

import argparse
import re
import unicodedata
from pathlib import Path

OUTPUT_DIR = Path("output/fasta")

RANK_PREFIXES = {
    "k": "kingdom",
    "p": "phylum",
    "c": "class",
    "o": "order",
    "f": "family",
    "g": "genus",
    "s": "species",
}

RANK_ORDER = ["kingdom", "phylum", "class", "order", "family", "genus", "species"]



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


def parse_header(line: str):
    """
    Parse '>ID;tax=k:V,p:V,...,s:V;' into (seq_id, ranks_dict).
    Returns (None, {}) if the tax= block is missing.
    """
    raw = line.lstrip(">").strip()
    parts = raw.split(";", 1)
    seq_id = parts[0].strip()

    ranks = {}
    tax_match = re.search(r"tax=([^;]+)", raw)
    if not tax_match:
        return seq_id, ranks

    for token in tax_match.group(1).split(","):
        token = token.strip()
        if ":" not in token:
            continue
        prefix, name = token.split(":", 1)
        prefix = prefix.strip()
        name = name.strip().replace("_", " ")
        rank = RANK_PREFIXES.get(prefix)
        if rank and name:
            ranks[rank] = name

    return seq_id, ranks


def get_scientific_name_and_rank(ranks: dict) -> tuple:
    for rank in reversed(RANK_ORDER):
        if rank in ranks:
            return sanitize(ranks[rank]), rank
    return "", ""


def get_higher_classification(ranks: dict) -> str:
    return ";".join(sanitize(ranks[r]) for r in RANK_ORDER if r in ranks)


def build_header(seq_id: str, ranks: dict, dataset: str, target_gene: str) -> str:
    scientific_name, taxon_rank = get_scientific_name_and_rank(ranks)
    higher_classification       = get_higher_classification(ranks)

    def r(rank: str) -> str:
        return sanitize(ranks.get(rank, ""))

    fields = [
        seq_id,
        seq_id,          # accessionNumber (no separate accession in source)
        scientific_name,
        "",              # decimalLatitude
        "",              # decimalLongitude
        "",              # typeStatus
        "",              # catalogueNumber
        "",              # identifiedBy
        taxon_rank,
        "",              # country
        "",              # locality
        "",              # basisOfRecord
        higher_classification,
        dataset,
        target_gene,
        "",              # domain
        r("kingdom"),
        r("phylum"),
        r("class"),
        r("order"),
        r("family"),
        r("genus"),
        r("species"),
    ]
    return "|".join(strip_delimiters(f) for f in fields)


def main():
    parser = argparse.ArgumentParser(description="Convert ITS2 Global FASTA to normalised header format")
    parser.add_argument("fasta_file",    help="Path to ITS2 Global FASTA file")
    parser.add_argument("dataset",       help="Short dataset name for headers (e.g. its2_global)")
    parser.add_argument("--target-gene", required=True, help="Target gene label (e.g. ITS2)")
    parser.add_argument("--output-dir",  default=None,  help="Directory to write output FASTA (default: output/fasta)")
    args = parser.parse_args()

    OUTPUT_DIR  = Path(args.output_dir) if args.output_dir else Path("output/fasta")
    fasta_path  = Path(args.fasta_file)
    target_gene = args.target_gene
    output_path = OUTPUT_DIR / f"{args.dataset}.fasta"

    if not fasta_path.exists():
        raise FileNotFoundError(f"FASTA not found: {fasta_path}")

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    print(f"Converting {fasta_path.name} …")
    written = skipped = 0
    seq_id = None
    ranks = {}
    seq_lines = []

    def flush():
        nonlocal written, skipped
        if seq_id is None or not seq_lines:
            return
        if not ranks:
            skipped += 1
            return
        header = build_header(seq_id, ranks, args.dataset, target_gene)
        out.write(f">{header}\n{''.join(seq_lines).upper()}\n")
        written += 1

    with open(fasta_path, encoding="utf-8") as inp, \
         open(output_path, "w", encoding="utf-8") as out:

        for line in inp:
            line = line.rstrip("\n")
            if line.startswith(">"):
                flush()
                seq_lines = []
                seq_id, ranks = parse_header(line)
            else:
                seq_lines.append(line.strip())

        flush()

    print(f"Done — {written:,} sequences written to {output_path}")
    if skipped:
        print(f"  {skipped:,} skipped (no parseable taxonomy)")


if __name__ == "__main__":
    main()
