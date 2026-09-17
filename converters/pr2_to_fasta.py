"""
Convert a PR2 merged xlsx file to the normalised FASTA header format.

Reads pr2_version_*_merged.xlsx and writes a FASTA whose headers follow the
same pipe-separated field order as dwc_to_fasta.py, midori2_to_fasta.py, and
mitofish_to_fasta.py.

Header format:
  ID | accessionNumber | scientificName | decimalLatitude | decimalLongitude |
  typeStatus | catalogueNumber | identifiedBy | taxonRank | country | locality |
  basisOfRecord | higherClassification | dataset | targetGene

Usage:
    python3 converters/pr2_to_fasta.py <xlsx_file> <dataset_shortname> --target-gene 18s

Example:
    python3 converters/pr2_to_fasta.py \\
        source-data/PR2/pr2_version_5.1.0_merged.xlsx \\
        pr2_18s --target-gene 18s
"""

import argparse
import re
import unicodedata
from pathlib import Path

import pandas as pd

OUTPUT_DIR = Path("output/fasta")

TAXON_FIELDS = ["domain", "supergroup", "division", "subdivision",
                "class", "order", "family", "genus", "species"]

NEEDED_COLS = [
    "pr2_main_id", "genbank_accession", "gene", "sequence",
    "pr2_latitude", "pr2_longitude", "gb_note", "gb_country", "pr2_sample_type",
] + TAXON_FIELDS

# Map user target-gene label to PR2 gene column value
GENE_NAME_MAP = {
    "18s": "18S_rRNA",
    "ssu_rrna_18s_eukaryotic": "18S_rRNA",
    "16s": "16S_rRNA",
}



# '|' separates header fields and '>' starts a FASTA record, so neither may appear
# inside a field value. Source data does contain them: NBDL's identifiedBy lists
# multiple collectors as "Pogonoski | Russell", which silently shifted every later
# field for 28 records until a header-width check caught it.
def strip_delimiters(value: str) -> str:
    return str(value).replace("|", "/").replace(">", "") if value else ""

def normalise_gene(target_gene: str) -> str:
    return GENE_NAME_MAP.get(target_gene.lower(), target_gene)


def sanitize(value) -> str:
    if not value or (isinstance(value, float)):
        return ""
    s = str(value)
    normalized = unicodedata.normalize("NFD", s)
    ascii_only = normalized.encode("ascii", "ignore").decode("ascii")
    ascii_only = strip_delimiters(ascii_only)
    return re.sub(r"\s+", "_", ascii_only).strip("_")


def clean_taxon(value) -> str:
    """Strip organelle suffixes (:nucl, :plas) from a taxonomy field."""
    if not value or (isinstance(value, float)):
        return ""
    return str(value).replace(":nucl", "").replace(":plas", "")


def get_taxon_rank_and_name(row) -> tuple[str, str]:
    """
    Find the lowest-resolved rank (no _X placeholder).
    A species ending in _sp. where genus matches is treated as genus-level.
    """
    for field in reversed(TAXON_FIELDS):
        val = row.get(field)
        if not val or (isinstance(val, float)):
            continue
        cleaned = clean_taxon(val)
        if "_X" in cleaned:
            continue
        if field == "species" and cleaned.endswith("_sp."):
            genus = clean_taxon(row.get("genus", ""))
            if genus and cleaned.replace("_sp.", "") == genus:
                return "genus", genus
        return field, cleaned
    return "", ""


def get_higher_classification(row) -> str:
    """Semicolon-separated taxonomy from domain down, skipping empty/_X fields."""
    parts = []
    for field in TAXON_FIELDS:
        val = clean_taxon(row.get(field, ""))
        if val and "_X" not in val:
            parts.append(val.replace(" ", "_"))
    return ";".join(parts)


def get_type_status(gb_note) -> str:
    if not gb_note or (isinstance(gb_note, float)):
        return ""
    note = str(gb_note).lower()
    for status in ("type strain", "holotype", "isotype", "lectotype", "neotype", "epitype"):
        if note.startswith(status) if status == "type strain" else status in note:
            return status
    return ""


def get_basis_of_record(pr2_sample_type) -> str:
    if not pr2_sample_type or (isinstance(pr2_sample_type, float)):
        return ""
    val = str(pr2_sample_type).lower()
    if val in ("isolate", "culture"):
        return "specimen"
    return "materialSample"


def split_country(gb_country) -> tuple[str, str]:
    if not gb_country or (isinstance(gb_country, float)):
        return "", ""
    parts = str(gb_country).split(":", 1)
    country  = sanitize(parts[0]) if len(parts) > 0 else ""
    locality = sanitize(parts[1]) if len(parts) > 1 else ""
    return country, locality


def get_rank_field(row, field: str) -> str:
    val = clean_taxon(row.get(field, ""))
    if not val or "_X" in val:
        return ""
    return sanitize(val)


def build_header(row, dataset: str, target_gene: str) -> str:
    taxon_rank, scientific_name = get_taxon_rank_and_name(row)
    country, locality           = split_country(row.get("gb_country"))

    fields = [
        str(row["pr2_main_id"]),
        sanitize(row.get("genbank_accession")),
        sanitize(scientific_name),
        sanitize(row.get("pr2_latitude")),
        sanitize(row.get("pr2_longitude")),
        get_type_status(row.get("gb_note")),
        "",                                         # catalogueNumber not in PR2
        "",                                         # identifiedBy not in PR2
        taxon_rank,
        country,
        locality,
        get_basis_of_record(row.get("pr2_sample_type")),
        get_higher_classification(row),
        dataset,
        target_gene,
        get_rank_field(row, "domain"),
        get_rank_field(row, "supergroup"),          # supergroup → kingdom
        get_rank_field(row, "division"),            # division → phylum
        get_rank_field(row, "class"),
        get_rank_field(row, "order"),
        get_rank_field(row, "family"),
        get_rank_field(row, "genus"),
        get_rank_field(row, "species"),
    ]
    return "|".join(strip_delimiters(f) for f in fields)


def main():
    parser = argparse.ArgumentParser(description="Convert PR2 xlsx to normalised FASTA format")
    parser.add_argument("xlsx_file",  help="Path to PR2 merged xlsx file")
    parser.add_argument("dataset",    help="Short dataset name for headers (e.g. pr2_18s)")
    parser.add_argument("--target-gene", required=True, help="Target gene to extract (e.g. 18s)")
    parser.add_argument("--output-dir",  default=None, help="Directory to write output FASTA (default: output/fasta)")
    args = parser.parse_args()

    OUTPUT_DIR  = Path(args.output_dir) if args.output_dir else Path("output/fasta")
    xlsx_path   = Path(args.xlsx_file)
    target_gene = args.target_gene
    gene_name   = normalise_gene(target_gene)
    output_path = OUTPUT_DIR / f"{args.dataset}.fasta"

    if not xlsx_path.exists():
        raise FileNotFoundError(f"xlsx not found: {xlsx_path}")

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    print(f"Reading {xlsx_path.name} (columns: {len(NEEDED_COLS)}) …")
    df = pd.read_excel(xlsx_path, usecols=NEEDED_COLS, dtype=str, engine="openpyxl")
    print(f"  {len(df):,} rows loaded")

    df = df[df["gene"] == gene_name]
    print(f"  {len(df):,} rows after filtering for {gene_name}")

    print(f"Writing {output_path} …")
    written = skipped = 0
    with open(output_path, "w", encoding="utf-8") as out:
        for _, row in df.iterrows():
            seq = str(row.get("sequence", "")).strip().upper()
            if not seq or seq == "NAN":
                skipped += 1
                continue
            header = build_header(row, args.dataset, target_gene)
            out.write(f">{header}\n{seq}\n")
            written += 1

    print(f"\nDone — {written:,} sequences written to {output_path}")
    if skipped:
        print(f"  {skipped} skipped (empty sequence)")


if __name__ == "__main__":
    main()
