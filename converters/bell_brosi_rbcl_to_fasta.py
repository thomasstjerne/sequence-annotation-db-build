"""
Convert the Bell & Brosi rbcL reference library FASTA to the normalised header format.

Source FASTA header format (DADA2 taxonomy string, no separate accession):
  >k__Kingdom_taxid;p__Phylum_taxid;c__sub__ClassName_taxid;o__Order_taxid;f__Family_taxid;g__Genus_taxid;s__Species name_taxid;

Taxids are numeric suffixes appended with an underscore to each rank name.
Class nodes have an optional 'sub__' infix: 'c__sub__asterids_71274'.

Output header format (pipe-separated, same field order as all other conversion scripts):
  ID | accessionNumber | scientificName | decimalLatitude | decimalLongitude |
  typeStatus | catalogueNumber | identifiedBy | taxonRank | country | locality |
  basisOfRecord | higherClassification | dataset | targetGene |
  domain | kingdom | phylum | class | order | family | genus | species

IDs are assigned as sequential integers (1, 2, 3, …).

Usage:
    python3 converters/bell_brosi_rbcl_to_fasta.py <fasta_file> <dataset_shortname> --target-gene rbcL

Example:
    python3 converters/bell_brosi_rbcl_to_fasta.py \\
        source-data/bell_brosi_rbcl/28757943 \\
        bell_brosi_rbcl --target-gene rbcL
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

# Matches trailing _<digits> taxid appended to every rank name
TAXID_SUFFIX = re.compile(r"_\d+$")



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


def strip_taxid(name: str) -> str:
    return TAXID_SUFFIX.sub("", name)


def parse_taxonomy(header: str) -> dict:
    """Parse 'k__X_tid;p__Y_tid;...;' into {rank: name}, skipping empty nodes."""
    ranks = {}
    for node in header.rstrip(";").split(";"):
        node = node.strip()
        if not node:
            continue
        for prefix, rank in RANK_PREFIXES.items():
            if node.startswith(prefix):
                raw = node[len(prefix):]
                # Handle optional 'sub__' infix on class nodes (c__sub__asterids_71274)
                if raw.startswith("sub__"):
                    raw = raw[len("sub__"):]
                name = strip_taxid(raw).strip()
                if name:
                    ranks[rank] = name
                break
    return ranks


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
    parser = argparse.ArgumentParser(description="Convert Bell & Brosi rbcL FASTA to normalised header format")
    parser.add_argument("fasta_file",    help="Path to rbcL DADA2 FASTA file")
    parser.add_argument("dataset",       help="Short dataset name for headers (e.g. bell_brosi_rbcl)")
    parser.add_argument("--target-gene", required=True, help="Target gene label (e.g. rbcL)")
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
    tax_header = ""
    seq_lines = []
    counter = 0

    def flush():
        nonlocal written, skipped
        if seq_id is None or not seq_lines:
            return
        ranks = parse_taxonomy(tax_header)
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
                counter += 1
                seq_id = str(counter)
                tax_header = line[1:].strip()
            else:
                seq_lines.append(line.strip())

        flush()

    print(f"Done — {written:,} sequences written to {output_path}")
    if skipped:
        print(f"  {skipped:,} skipped (no parseable taxonomy)")


if __name__ == "__main__":
    main()
