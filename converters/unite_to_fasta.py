"""
Convert a UNITE general release FASTA to the normalised header format.

UNITE FASTA header format:
  >scientificName|accession|SH_ID|type|k__kingdom;p__phylum;c__class;o__order;f__family;g__genus;s__species

Output header format (pipe-separated, same field order as all other conversion scripts):
  ID | accessionNumber | scientificName | decimalLatitude | decimalLongitude |
  typeStatus | catalogueNumber | identifiedBy | taxonRank | country | locality |
  basisOfRecord | higherClassification | dataset | targetGene

The SH (Species Hypothesis) ID is used as the record ID.
Incertae sedis placeholders and _sp suffixes are treated as unresolved.

Usage:
    python3 converters/unite_to_fasta.py <fasta_file> <dataset_shortname> --target-gene its

Example:
    python3 converters/unite_to_fasta.py \\
        source-data/UNITE/sh_general_release_dynamic_s_all_19.02.2025.fasta \\
        unite_its --target-gene its
"""

import argparse
import re
import unicodedata
from pathlib import Path

OUTPUT_DIR = Path("output/fasta")

RANK_PREFIXES = {
    "k__": "kingdom",
    "p__": "phylum",
    "c__": "class",
    "o__": "order",
    "f__": "family",
    "g__": "genus",
    "s__": "species",
}

RANK_ORDER = ["kingdom", "phylum", "class", "order", "family", "genus", "species"]

UNRESOLVED = re.compile(r"_Incertae_sedis|_kgd_|_phy_|_cls_|_ord_|_fam_|_gen_")



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


def is_unresolved(name: str) -> bool:
    return bool(UNRESOLVED.search(name)) or name.endswith("_sp")


def parse_taxonomy(tax_str: str) -> dict[str, str]:
    """Parse 'k__X;p__Y;...' into {rank: name}, skipping unresolved nodes."""
    ranks = {}
    for node in tax_str.split(";"):
        node = node.strip()
        for prefix, rank in RANK_PREFIXES.items():
            if node.startswith(prefix):
                name = node[len(prefix):]
                if name and not is_unresolved(name):
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


def parse_header(line: str) -> tuple[str, str, str, str, str, str, dict]:
    """
    Returns: sh_id, accession, scientific_name, taxon_rank, higher_classification, sh_type, ranks
    """
    raw = line.lstrip(">").strip()
    parts = raw.split("|")

    accession    = parts[1].strip() if len(parts) > 1 else ""
    sh_id        = parts[2].strip() if len(parts) > 2 else raw
    sh_type      = parts[3].strip() if len(parts) > 3 else ""
    tax_str      = parts[4].strip() if len(parts) > 4 else ""

    ranks                         = parse_taxonomy(tax_str)
    scientific_name, taxon_rank   = get_scientific_name_and_rank(ranks)
    higher_classification         = get_higher_classification(ranks)

    return sh_id, accession, scientific_name, taxon_rank, higher_classification, sh_type, ranks


def build_header(sh_id, accession, scientific_name, taxon_rank,
                 higher_classification, ranks, dataset, target_gene) -> str:
    def r(rank: str) -> str:
        return sanitize(ranks.get(rank, ""))

    fields = [
        sh_id,
        accession,
        scientific_name,
        "",    # decimalLatitude
        "",    # decimalLongitude
        "",    # typeStatus
        "",    # catalogueNumber
        "",    # identifiedBy
        taxon_rank,
        "",    # country
        "",    # locality
        "",    # basisOfRecord
        higher_classification,
        dataset,
        target_gene,
        "",           # domain
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
    parser = argparse.ArgumentParser(description="Convert UNITE FASTA to normalised header format")
    parser.add_argument("fasta_file",  help="Path to extracted UNITE FASTA file")
    parser.add_argument("dataset",     help="Short dataset name for headers (e.g. unite_its)")
    parser.add_argument("--target-gene", required=True, help="Target gene label (e.g. its)")
    parser.add_argument("--output-dir",  default=None, help="Directory to write output FASTA (default: output/fasta)")
    args = parser.parse_args()

    OUTPUT_DIR  = Path(args.output_dir) if args.output_dir else Path("output/fasta")
    fasta_path  = Path(args.fasta_file)
    target_gene = args.target_gene.lower()
    output_path = OUTPUT_DIR / f"{args.dataset}.fasta"

    if not fasta_path.exists():
        raise FileNotFoundError(f"FASTA not found: {fasta_path}")

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    print(f"Converting {fasta_path.name} …")
    written = 0

    with open(fasta_path, encoding="utf-8") as inp, \
         open(output_path, "w", encoding="utf-8") as out:

        sh_id = accession = scientific_name = taxon_rank = higher_classification = ""
        cur_ranks: dict = {}
        seq_lines = []

        def flush():
            nonlocal written
            if sh_id and seq_lines:
                header = build_header(sh_id, accession, scientific_name, taxon_rank,
                                      higher_classification, cur_ranks, args.dataset, target_gene)
                out.write(f">{header}\n{''.join(seq_lines).upper()}\n")
                written += 1

        for line in inp:
            line = line.rstrip("\n")
            if line.startswith(">"):
                flush()
                seq_lines = []
                sh_id, accession, scientific_name, taxon_rank, higher_classification, _, cur_ranks = \
                    parse_header(line)
            else:
                seq_lines.append(line.strip())

        flush()

    print(f"Done — {written:,} sequences written to {output_path}")


if __name__ == "__main__":
    main()
