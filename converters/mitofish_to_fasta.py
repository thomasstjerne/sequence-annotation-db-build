"""
Convert MitoFish to the normalised FASTA header format.

MitoFish consists of a gzipped FASTA (bare accession headers) plus a set of
Parquet tables that map accessions to gene annotations and taxonomy.

Only sequences with an annotation for the requested target gene are included.
Where the annotated region is on the minus strand, the subsequence is
reverse-complemented before output.

Output header format (pipe-separated, same field order as dwc_to_fasta.py):
  ID | accessionNumber | scientificName | decimalLatitude | decimalLongitude |
  typeStatus | catalogueNumber | identifiedBy | taxonRank | country | locality |
  basisOfRecord | higherClassification | dataset | targetGene

Fields absent from MitoFish (lat/lon, typeStatus, etc.) are left empty.
higherClassification uses the ranks available in MitoFish: class, order,
family, genus, species (kingdom/phylum are not stored as MitoFish is fish-only).

Usage:
    python3 converters/mitofish_to_fasta.py <fasta_gz> <tables_dir> <dataset_shortname> --target-gene 12s

Example:
    python3 converters/mitofish_to_fasta.py \\
        source-data/mitofish/mitofishdb.fa.gz \\
        source-data/mitofish/tables \\
        mitofish_12s --target-gene 12s
"""

import argparse
import gzip
import re
from pathlib import Path

import duckdb

OUTPUT_DIR = Path("output/fasta")

COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")

# Map user-supplied target gene label to the gene name used in seq_annotation
GENE_NAME_MAP = {
    "12s": "12S_rRNA",
    "ssu12smitochondrial": "12S_rRNA",
    "16s": "16S_rRNA",
    "coi": "COI",
}



# '|' separates header fields and '>' starts a FASTA record, so neither may appear
# inside a field value. Source data does contain them: NBDL's identifiedBy lists
# multiple collectors as "Pogonoski | Russell", which silently shifted every later
# field for 28 records until a header-width check caught it.
def strip_delimiters(value: str) -> str:
    return str(value).replace("|", "/").replace(">", "") if value else ""

def revcomp(seq: str) -> str:
    return seq.translate(COMPLEMENT)[::-1]


def normalise_gene(target_gene: str) -> str:
    key = target_gene.lower().replace("_rrna", "").replace("-", "").replace("_", "")
    return GENE_NAME_MAP.get(key, target_gene)


def load_annotations(tables_dir: Path, gene_name: str) -> dict[str, tuple[int, int, str]]:
    """Return {accession: (start, end, strand)} for all accessions annotated with gene_name."""
    path = tables_dir / "seq_annotation.parquet"
    con = duckdb.connect()
    rows = con.execute(
        f'SELECT accession, "start", "end", strand FROM \'{path}\' WHERE gene = ?',
        [gene_name],
    ).fetchall()
    return {row[0]: (row[1], row[2], row[3]) for row in rows}


def load_taxonomy(tables_dir: Path, accessions: set[str]) -> dict[str, dict]:
    """
    Return {accession: {class, order, family, genus, species,
                        scientificName, taxonRank, higherClassification}}
    for the given set of accessions.
    """
    con = duckdb.connect()

    seq_tax   = tables_dir / "seq_taxonid.parquet"
    tax_sp    = tables_dir / "taxonid_speciesid.parquet"
    sp_lin    = tables_dir / "speciesid_lineageid.parquet"
    lin_name  = tables_dir / "taxonid_name.parquet"

    rows = con.execute(f"""
        SELECT
            s.accession,
            n.lineage_rank,
            n.lineage_name
        FROM '{seq_tax}' s
        JOIN '{tax_sp}'  ts  ON s.taxon_id   = ts.taxon_id
        JOIN '{sp_lin}'  sl  ON ts.species_id = sl.species_id
        JOIN '{lin_name}' n  ON sl.lineage_id = n.lineage_id
        WHERE n.lineage_rank IN ('class','order','family','genus','species','subspecies')
        ORDER BY s.accession, sl.rank_order
    """).fetchall()

    STANDARD_RANKS = ["class", "order", "family", "genus", "species"]

    taxonomy: dict[str, dict] = {}
    for accession, rank, name in rows:
        if accession not in taxonomy:
            taxonomy[accession] = {}
        taxonomy[accession][rank] = name

    result = {}
    for accession, ranks in taxonomy.items():
        # Build higherClassification from standard ranks in order
        parts = [ranks[r] for r in STANDARD_RANKS if r in ranks]
        # Prefer species name; fall back to lowest available rank
        sci_name = ranks.get("species") or ranks.get("genus") or ""
        taxon_rank = next(
            (r for r in reversed(STANDARD_RANKS) if r in ranks), ""
        )
        # Replace spaces with underscores for header compatibility
        hc = ";".join(p.replace(" ", "_") for p in parts)
        result[accession] = {
            "scientificName":       sci_name.replace(" ", "_"),
            "taxonRank":            taxon_rank,
            "higherClassification": hc,
            "class":                ranks.get("class", "").replace(" ", "_"),
            "order":                ranks.get("order", "").replace(" ", "_"),
            "family":               ranks.get("family", "").replace(" ", "_"),
            "genus":                ranks.get("genus", "").replace(" ", "_"),
            "species":              ranks.get("species", "").replace(" ", "_"),
        }

    return result


def build_header(accession: str, tax: dict, dataset: str, target_gene: str) -> str:
    fields = [
        accession,
        accession,                              # accessionNumber same as ID
        tax.get("scientificName", ""),
        "",                                     # decimalLatitude
        "",                                     # decimalLongitude
        "",                                     # typeStatus
        "",                                     # catalogueNumber
        "",                                     # identifiedBy
        tax.get("taxonRank", ""),
        "",                                     # country
        "",                                     # locality
        "",                                     # basisOfRecord
        tax.get("higherClassification", ""),
        dataset,
        target_gene,
        "Eukaryota",                            # domain
        "Animalia",                             # kingdom
        "Chordata",                             # phylum
        tax.get("class", ""),
        tax.get("order", ""),
        tax.get("family", ""),
        tax.get("genus", ""),
        tax.get("species", ""),
    ]
    return "|".join(strip_delimiters(f) for f in fields)


def main():
    parser = argparse.ArgumentParser(description="Convert MitoFish to normalised FASTA format")
    parser.add_argument("fasta_gz",    help="Path to mitofishdb.fa.gz")
    parser.add_argument("tables_dir",  help="Directory containing the extracted MitoFish parquet tables")
    parser.add_argument("dataset",     help="Short dataset name for headers (e.g. mitofish_12s)")
    parser.add_argument("--target-gene", required=True, help="Target gene to extract (e.g. 12s)")
    parser.add_argument("--output-dir",  default=None, help="Directory to write output FASTA (default: output/fasta)")
    args = parser.parse_args()

    OUTPUT_DIR  = Path(args.output_dir) if args.output_dir else Path("output/fasta")
    fasta_path  = Path(args.fasta_gz)
    tables_dir  = Path(args.tables_dir)
    target_gene = args.target_gene
    gene_name   = normalise_gene(target_gene)
    output_path = OUTPUT_DIR / f"{args.dataset}.fasta"

    if not fasta_path.exists():
        raise FileNotFoundError(f"FASTA not found: {fasta_path}")
    if not tables_dir.exists():
        raise FileNotFoundError(f"Tables directory not found: {tables_dir}")

    OUTPUT_DIR.mkdir(exist_ok=True)

    print(f"Loading {gene_name} annotations …")
    annotations = load_annotations(tables_dir, gene_name)
    print(f"  {len(annotations):,} annotated accessions")

    print("Loading taxonomy …")
    taxonomy = load_taxonomy(tables_dir, set(annotations.keys()))
    print(f"  {len(taxonomy):,} taxonomy records")

    print(f"Streaming {fasta_path.name} …")
    written = skipped_no_annotation = skipped_no_sequence = 0

    with gzip.open(fasta_path, "rt", encoding="utf-8") as inp, \
         open(output_path, "w", encoding="utf-8") as out:

        current_id = None
        seq_lines  = []

        def flush():
            nonlocal written, skipped_no_annotation, skipped_no_sequence
            if not current_id:
                return
            if current_id not in annotations:
                skipped_no_annotation += 1
                return
            seq = "".join(seq_lines).upper()
            if not seq:
                skipped_no_sequence += 1
                return

            start, end, strand = annotations[current_id]
            # GenBank coordinates are 1-based inclusive
            subseq = seq[start - 1:end]
            if strand == "-1":
                subseq = revcomp(subseq)

            tax    = taxonomy.get(current_id, {})
            header = build_header(current_id, tax, args.dataset, target_gene)
            out.write(f">{header}\n{subseq}\n")
            written += 1

        for line in inp:
            line = line.rstrip("\n")
            if line.startswith(">"):
                flush()
                current_id = line[1:].split()[0]
                seq_lines  = []
            else:
                seq_lines.append(line.strip())

        flush()  # last record

    print(f"\nDone — {written:,} sequences written to {output_path}")
    print(f"  {skipped_no_annotation:,} skipped (no {gene_name} annotation)")
    if skipped_no_sequence:
        print(f"  {skipped_no_sequence:,} skipped (empty sequence)")


if __name__ == "__main__":
    main()
