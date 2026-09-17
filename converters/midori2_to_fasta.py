"""
Convert a MIDORI2 BLAST/UNIQ FASTA to the normalised header format.

MIDORI2 header format:
  >accession.version.start.end###root_1;Taxon_taxid;...;Genus_species_taxid

Output header format (pipe-separated, same field order as dwc_to_fasta.py):
  ID | accessionNumber | scientificName | decimalLatitude | decimalLongitude |
  typeStatus | catalogueNumber | identifiedBy | taxonRank | country | locality |
  basisOfRecord | higherClassification | dataset | targetGene |
  domain | kingdom | phylum | class | order | family | genus | species

Fields absent from MIDORI2 (lat/lon, typeStatus, etc.) are left empty.

Usage:
    python3 converters/midori2_to_fasta.py <fasta_file> <dataset_shortname> --target-gene 12s

Example:
    python3 converters/midori2_to_fasta.py \\
        source-data/midori2/MIDORI2_UNIQ_NUC_GB269_srRNA_BLAST/MIDORI2_UNIQ_NUC_GB269_srRNA_BLAST.fasta \\
        midori2_12s --target-gene 12s
"""

import argparse
import re
from pathlib import Path

OUTPUT_DIR = Path("output/fasta")

RANK_PREFIXES = {"kingdom", "phylum", "class", "order", "family", "genus", "suborder",
                 "subclass", "subfamily", "superfamily", "tribe", "subtribe", "infraorder"}



# '|' separates header fields and '>' starts a FASTA record, so neither may appear
# inside a field value. Source data does contain them: NBDL's identifiedBy lists
# multiple collectors as "Pogonoski | Russell", which silently shifted every later
# field for 28 records until a header-width check caught it.
def strip_delimiters(value: str) -> str:
    return str(value).replace("|", "/").replace(">", "") if value else ""

def parse_node(node: str) -> tuple[str, str, str]:
    """
    Parse a single taxonomy node into (rank, name, taxid).

    Examples:
      'Eukaryota_2759'              → ('', 'Eukaryota', '2759')
      'order_Carcharhiniformes_3639' → ('order', 'Carcharhiniformes', '3639')
      'Paravannella_minima_1443144'  → ('', 'Paravannella minima', '1443144')
    """
    parts = node.split("_")

    # Strip trailing numeric taxid
    if parts and parts[-1].isdigit():
        taxid = parts.pop()
    else:
        taxid = ""

    # Strip leading rank prefix
    rank = ""
    if parts and parts[0].lower() in RANK_PREFIXES:
        rank = parts.pop(0).lower()

    name = " ".join(parts)
    return rank, name, taxid


def infer_rank(rank_prefix: str, name: str) -> str:
    """Infer taxon rank from the rank prefix or the number of words in the name."""
    if rank_prefix:
        return rank_prefix
    if len(name.split()) == 2:
        return "species"
    return ""


def parse_taxonomy(taxonomy_str: str) -> tuple[str, str, str, str, str, str]:
    """
    Parse the MIDORI2 taxonomy string (after ###) into:
      - scientificName, taxonRank, higherClassification
      - family, genus, species  (extracted from node structure)
    """
    nodes = [n for n in taxonomy_str.split(";") if n and not n.startswith("root_")]
    if not nodes:
        return "", "", "", "", "", ""

    parsed = [(parse_node(n)[0], parse_node(n)[1]) for n in nodes]  # (rank_prefix, name)
    names = [name for _, name in parsed if name]

    if not names:
        return "", "", "", "", "", ""

    scientific_name = names[-1]
    last_rank_prefix = parsed[-1][0]
    taxon_rank = infer_rank(last_rank_prefix, scientific_name)
    higher_classification = ";".join(names)

    # Extract family: node with explicit rank "family", or name ending in -idae/-aceae
    family = ""
    genus = ""
    species = ""
    for rank_prefix, name in parsed:
        if rank_prefix == "family":
            family = name
        elif not rank_prefix and (name.endswith("idae") or name.endswith("aceae")):
            family = name
        elif rank_prefix == "genus":
            genus = name
        elif rank_prefix == "subgenus":
            pass
    # Genus: second-to-last node if the last is a species binomial
    if not genus and len(names) >= 2:
        candidate = names[-2]
        if len(names[-1].split()) == 2:
            genus = candidate
    # Species: last node only if it's a binomial
    if len(scientific_name.split()) == 2:
        species = scientific_name

    return scientific_name, taxon_rank, higher_classification, family, genus, species


def parse_header(header_line: str) -> tuple[str, str, str, str, str, str, str, str]:
    """
    Parse a MIDORI2 FASTA header line into:
      seq_id, accession_number, scientific_name, taxon_rank, higher_classification,
      family, genus, species
    """
    raw = header_line.lstrip(">").strip()

    if "###" in raw:
        seq_id, taxonomy_str = raw.split("###", 1)
    else:
        seq_id = raw
        taxonomy_str = ""

    # GenBank accession is the first two dot-separated components (e.g. MH910097.1)
    acc_parts = seq_id.split(".")
    accession_number = ".".join(acc_parts[:2]) if len(acc_parts) >= 2 else seq_id

    scientific_name, taxon_rank, higher_classification, family, genus, species = \
        parse_taxonomy(taxonomy_str)

    return seq_id, accession_number, scientific_name, taxon_rank, higher_classification, \
        family, genus, species


def build_header(seq_id, accession_number, scientific_name, taxon_rank,
                 higher_classification, family, genus, species,
                 dataset, target_gene) -> str:
    def s(v):
        return v.replace(" ", "_") if v else ""

    fields = [
        seq_id,
        accession_number,
        s(scientific_name),
        "",   # decimalLatitude
        "",   # decimalLongitude
        "",   # typeStatus
        "",   # catalogueNumber
        "",   # identifiedBy
        taxon_rank,
        "",   # country
        "",   # locality
        "",   # basisOfRecord
        s(higher_classification),
        dataset,
        target_gene,
        "Eukaryota",    # domain
        "",             # kingdom
        "",             # phylum
        "",             # class
        "",             # order
        s(family),
        s(genus),
        s(species),
    ]
    return "|".join(strip_delimiters(f) for f in fields)


def main():
    parser = argparse.ArgumentParser(description="Convert MIDORI2 FASTA to normalised header format")
    parser.add_argument("fasta_file", help="Path to MIDORI2 FASTA file")
    parser.add_argument("dataset",    help="Short dataset name for headers (e.g. midori2_12s)")
    parser.add_argument("--target-gene", required=True, help="Target gene label (e.g. 12s)")
    parser.add_argument("--output-dir",  default=None, help="Directory to write output FASTA (default: output/fasta)")
    args = parser.parse_args()

    OUTPUT_DIR  = Path(args.output_dir) if args.output_dir else Path("output/fasta")
    fasta_path  = Path(args.fasta_file)
    target_gene = args.target_gene.lower()
    output_path = OUTPUT_DIR / f"{args.dataset}.fasta"

    if not fasta_path.exists():
        raise FileNotFoundError(f"FASTA file not found: {fasta_path}")

    OUTPUT_DIR.mkdir(exist_ok=True)

    print(f"Converting {fasta_path.name} …")
    written = 0

    with open(fasta_path, encoding="utf-8") as inp, \
         open(output_path, "w", encoding="utf-8") as out:

        seq_id = accession_number = scientific_name = taxon_rank = higher_classification = ""
        family = genus = species = ""
        sequence_lines = []

        def flush():
            nonlocal written
            if seq_id and sequence_lines:
                header = build_header(seq_id, accession_number, scientific_name,
                                      taxon_rank, higher_classification,
                                      family, genus, species,
                                      args.dataset, target_gene)
                out.write(f">{header}\n{''.join(sequence_lines).upper()}\n")
                written += 1

        for line in inp:
            line = line.rstrip("\n")
            if line.startswith(">"):
                flush()
                sequence_lines = []
                seq_id, accession_number, scientific_name, taxon_rank, higher_classification, \
                    family, genus, species = parse_header(line)
            else:
                sequence_lines.append(line.strip())

        flush()  # last record

    print(f"Done — {written:,} sequences written to {output_path}")


if __name__ == "__main__":
    main()
