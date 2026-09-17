"""
Convert a Darwin Core archive to a normalised FASTA file.

Reads dna.tsv/dna.txt and occurrence.tsv/occurrence.txt from an extracted DwC
archive and writes a FASTA file whose headers contain the fields defined in
tsvHeaders.js plus dataset shortname and target gene.

Both named-header archives (e.g. NBDL) and positional/headerless archives
(e.g. refSeq ITS, following the standard DwC meta.xml field indices) are
supported. Header presence is auto-detected from the first row.

Header format (pipe-separated):
  ID | accessionNumber | scientificName | decimalLatitude | decimalLongitude |
  typeStatus | catalogueNumber | identifiedBy | taxonRank | country | locality |
  basisOfRecord | higherClassification | dataset | targetGene |
  domain | kingdom | phylum | class | order | family | genus | species |
  domain | kingdom | phylum | class | order | family | genus | species

Usage:
    python3 converters/dwc_to_fasta.py <dwc_dir> <dataset_shortname> --target-gene 12s

Examples:
    python3 converters/dwc_to_fasta.py source-data/nbdl/extracted nbdl_12s --target-gene 12s
    python3 converters/dwc_to_fasta.py source-data/refSeq_its refseq_its --target-gene its
"""

import argparse
import csv
import re
import unicodedata
from pathlib import Path

OUTPUT_DIR = Path("output/fasta")

DNA_ID       = 0
DNA_GENE     = 1
DNA_SEQUENCE = 2

# Standard DwC column names for archives without a header row,
# ordered by meta.xml field index (0 = core id).
DWC_STANDARD_HEADERS = [
    "id", "occurrenceID", "scientificName", "higherClassification",
    "typeStatus", "decimalLatitude", "decimalLongitude", "country", "locality",
    "eventDate", "recordedBy", "identifiedBy", "catalogNumber",
    "associatedReferences", "basisOfRecord", "occurrenceRemarks",
    "taxonID", "taxonConceptID", "associatedSequences", "otherCatalogNumbers",
    "kingdom", "phylum", "class", "order", "family", "genus", "institutionCode",
]

# First-column values that indicate a header row is present
HEADER_SENTINEL = {"occurrenceId", "occurrenceID", "id", "ID", "coreid", "coreId"}



# '|' separates header fields and '>' starts a FASTA record, so neither may appear
# inside a field value. Source data does contain them: NBDL's identifiedBy lists
# multiple collectors as "Pogonoski | Russell", which silently shifted every later
# field for 28 records until a header-width check caught it.
def strip_delimiters(value: str) -> str:
    return str(value).replace("|", "/").replace(">", "") if value else ""

def sanitize(value: str) -> str:
    """ASCII-safe, whitespace→underscore, for BLAST/vsearch header compatibility."""
    if not value:
        return ""
    normalized = unicodedata.normalize("NFD", value)
    ascii_only = normalized.encode("ascii", "ignore").decode("ascii")
    ascii_only = strip_delimiters(ascii_only)
    return re.sub(r"\s+", "_", ascii_only).strip("_")


def get_taxon_rank(row: dict) -> str:
    genus    = row.get("genus", "").strip()
    specific = row.get("specificEpiphet", "").strip()
    sci      = row.get("scientificName", "").strip()
    # Species: explicit specificEpiphet, or scientificName has ≥2 words starting with genus
    if genus and specific:
        return "species"
    if genus and sci and len(sci.split()) >= 2 and sci.startswith(genus):
        return "species"
    if genus:
        return "genus"
    if row.get("family", "").strip():
        return "family"
    if row.get("order", "").strip():
        return "order"
    if row.get("class", "").strip():
        return "class"
    if row.get("phylum", "").strip():
        return "phylum"
    if row.get("kingdom", "").strip():
        return "kingdom"
    return ""


def build_higher_classification(row: dict) -> str:
    """Semicolon-separated taxonomy from kingdom down to species binomial."""
    parts = []
    for rank in ("kingdom", "phylum", "class", "order", "family", "genus"):
        val = row.get(rank, "").strip()
        if val:
            parts.append(val)
    genus    = row.get("genus", "").strip()
    specific = row.get("specificEpiphet", "").strip()
    sci      = row.get("scientificName", "").strip()
    if genus and specific:
        parts.append(f"{genus}_{specific}")
    elif genus and sci and len(sci.split()) >= 2 and sci.startswith(genus):
        # Derive species binomial from scientificName when specificEpiphet absent
        parts.append("_".join(sci.split()[:2]))
    return ";".join(parts)


def get_species(row: dict) -> str:
    genus    = row.get("genus", "").strip()
    specific = row.get("specificEpiphet", "").strip()
    sci      = row.get("scientificName", "").strip()
    if genus and specific:
        return sanitize(f"{genus} {specific}")
    if genus and sci and len(sci.split()) >= 2 and sci.startswith(genus):
        return sanitize(" ".join(sci.split()[:2]))
    return ""


def load_sequences(dna_path: Path, filter_gene) -> dict[str, str]:
    sequences: dict[str, str] = {}
    with open(dna_path, newline="", encoding="utf-8") as f:
        reader = csv.reader(f, delimiter="\t")
        for row in reader:
            if len(row) < 3:
                continue
            if filter_gene is None or row[DNA_GENE].strip().lower() == filter_gene:
                sequences[row[DNA_ID].strip()] = row[DNA_SEQUENCE].strip().upper()
    return sequences


def load_occurrences(occ_path: Path) -> dict[str, dict]:
    occurrences: dict[str, dict] = {}
    with open(occ_path, newline="", encoding="utf-8") as f:
        first_col = f.readline().split("\t")[0].strip()
        f.seek(0)
        if first_col in HEADER_SENTINEL:
            reader = csv.DictReader(f, delimiter="\t")
        else:
            reader = csv.DictReader(f, fieldnames=DWC_STANDARD_HEADERS, delimiter="\t")
        for row in reader:
            # Support both named-header archives and standard DwC positional archives
            occ_id = (row.get("occurrenceId") or row.get("occurrenceID")
                      or row.get("id") or "").strip()
            if occ_id:
                occurrences[occ_id] = row
    return occurrences


def build_header(occ: dict, seq_id: str, dataset: str, target_gene: str) -> str:
    fields = [
        seq_id,
        sanitize(occ.get("otherCatalogNumbers", "")),
        sanitize(occ.get("scientificName", "")),
        sanitize(occ.get("decimalLatitude", "")),
        sanitize(occ.get("decimalLongitude", "")),
        sanitize(occ.get("typeStatus", "")),
        sanitize(occ.get("catalogNumber", "")),
        sanitize(occ.get("identifiedBy", "")),
        get_taxon_rank(occ),
        sanitize(occ.get("country", "")),
        sanitize(occ.get("locality", "")),
        sanitize(occ.get("basisOfRecord", "")),  # blank in NBDL
        sanitize(build_higher_classification(occ)),
        dataset,
        target_gene,
        "",                                       # domain
        sanitize(occ.get("kingdom", "")),
        sanitize(occ.get("phylum", "")),
        sanitize(occ.get("class", "")),
        sanitize(occ.get("order", "")),
        sanitize(occ.get("family", "")),
        sanitize(occ.get("genus", "")),
        get_species(occ),
    ]
    return "|".join(strip_delimiters(f) for f in fields)


def main():
    parser = argparse.ArgumentParser(description="Convert a DwC archive to a normalised FASTA")
    parser.add_argument("dwc_dir",  help="Path to extracted DwC archive directory")
    parser.add_argument("dataset",  help="Short dataset name for FASTA headers (e.g. nbdl_12s)")
    parser.add_argument("--target-gene", required=True, help="Vocabulary term for output FASTA headers (e.g. SSU_rRNA_12S_mitochondrial)")
    parser.add_argument("--filter-gene", default=None, help="Value to match against the targetGene field in the DWC source (e.g. 12s); if omitted, all sequences are included")
    parser.add_argument("--output-dir",  default=None, help="Directory to write output FASTA (default: output/fasta)")
    args = parser.parse_args()

    OUTPUT_DIR  = Path(args.output_dir) if args.output_dir else Path("output/fasta")
    target_gene = args.target_gene
    filter_gene = args.filter_gene.lower() if args.filter_gene else None
    dwc_dir     = Path(args.dwc_dir)
    output_path = OUTPUT_DIR / f"{args.dataset}.fasta"

    def find_file(stem: str) -> Path:
        for ext in (".tsv", ".txt"):
            p = dwc_dir / f"{stem}{ext}"
            if p.exists():
                return p
        raise FileNotFoundError(f"Could not find {stem}.tsv or {stem}.txt in {dwc_dir}")

    dna_path = find_file("dna")
    occ_path = find_file("occurrence")

    OUTPUT_DIR.mkdir(exist_ok=True)

    filter_desc = filter_gene if filter_gene else "all"
    print(f"Loading sequences (filter gene: {filter_desc}, output label: {target_gene}) …")
    sequences = load_sequences(dna_path, filter_gene)
    print(f"  {len(sequences):,} sequences found")

    print("Loading occurrences …")
    occurrences = load_occurrences(occ_path)
    print(f"  {len(occurrences):,} occurrence records")

    print(f"Writing {output_path} …")
    written = missing = 0
    with open(output_path, "w", encoding="utf-8") as out:
        for seq_id, sequence in sequences.items():
            occ = occurrences.get(seq_id)
            if not occ:
                print(f"  Warning: no occurrence record for {seq_id}")
                missing += 1
                continue
            header = build_header(occ, seq_id, args.dataset, target_gene)
            out.write(f">{header}\n{sequence}\n")
            written += 1

    print(f"\nDone — {written:,} sequences written to {output_path}")
    if missing:
        print(f"  {missing} sequences skipped (no occurrence record)")


if __name__ == "__main__":
    main()
