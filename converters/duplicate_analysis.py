"""
Analyse exact sequence duplicates in any normalised pipe-header FASTA.

Reports:
  - Total sequences vs unique sequences
  - How many sequences are exact duplicates
  - Whether duplicates are within-species or cross-species
  - Top duplicated sequences by copy count

With --write-deduped, writes a new FASTA that collapses within-species duplicates
(same sequence + same species → one entry) while keeping cross-species duplicates,
so that pickBestMatch can detect sequences shared across species and down-rank to genus.

Usage:
    python3 converters/duplicate_analysis.py --fasta output/fasta/its2_global.fasta
    python3 converters/duplicate_analysis.py --fasta output/fasta/pr2_18s.fasta --write-deduped output/fasta/pr2_18s_deduped.fasta
"""

import argparse
import sys
from pathlib import Path

import duckdb


HEADER_FIELDS = [
    'id', 'accessionNumber', 'scientificName',
    'decimalLatitude', 'decimalLongitude', 'typeStatus',
    'catalogueNumber', 'identifiedBy', 'taxonRank',
    'country', 'locality', 'basisOfRecord',
    'higherClassification', 'dataset', 'targetGene',
    'domain', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species',
]


def parse_fasta(path: Path):
    """Yield (header_fields_dict, raw_header_line, sequence) from a normalised pipe-header FASTA."""
    with open(path) as f:
        header = None
        raw_header = None
        bases = []
        for line in f:
            line = line.rstrip()
            if line.startswith('>'):
                if header is not None:
                    yield header, raw_header, ''.join(bases)
                parts = line[1:].split('|')
                header = {field: (parts[i] if i < len(parts) else '') for i, field in enumerate(HEADER_FIELDS)}
                raw_header = line
                bases = []
            else:
                bases.append(line)
        if header is not None:
            yield header, raw_header, ''.join(bases)


def main():
    parser = argparse.ArgumentParser(description="Analyse ITS2 FASTA for duplicate sequences")
    parser.add_argument('--fasta', default='output/fasta/its2_global.fasta',
                        help='Path to the ITS2 FASTA file')
    parser.add_argument('--write-deduped', metavar='PATH',
                        help='Write a within-species-deduplicated FASTA to this path')
    args = parser.parse_args()

    fasta_path = Path(args.fasta)
    if not fasta_path.exists():
        sys.exit(f"FASTA not found: {fasta_path}")

    print(f"Parsing {fasta_path} …", flush=True)
    rows = []
    # raw_records preserves (raw_header, sequence) for FASTA output
    raw_records = []
    for fields, raw_header, seq in parse_fasta(fasta_path):
        seq = seq.upper()
        rows.append({
            'id':             fields['id'],
            'species':        fields['species'],
            'genus':          fields['genus'],
            'family':         fields['family'],
            'scientificName': fields['scientificName'],
            'sequence':       seq,
        })
        raw_records.append((fields['species'], seq, raw_header))

    print(f"Loaded {len(rows):,} sequences — running analysis …", flush=True)

    import pandas as pd
    df = pd.DataFrame(rows)

    con = duckdb.connect()
    con.execute("CREATE TABLE seqs AS SELECT * FROM df")

    total = con.execute("SELECT COUNT(*) FROM seqs").fetchone()[0]
    unique_seqs = con.execute("SELECT COUNT(DISTINCT sequence) FROM seqs").fetchone()[0]
    duplicated = total - unique_seqs

    print(f"\n── Overview ──────────────────────────────────────────")
    print(f"  Total sequences:          {total:>10,}")
    print(f"  Unique sequences:         {unique_seqs:>10,}")
    print(f"  Duplicate copies:         {duplicated:>10,}  ({duplicated/total*100:.1f}% of total)")

    # sequences that appear more than once
    dup_groups = con.execute("""
        SELECT COUNT(*) AS n_groups, SUM(copies) AS n_copies
        FROM (
            SELECT sequence, COUNT(*) AS copies
            FROM seqs
            GROUP BY sequence
            HAVING COUNT(*) > 1
        )
    """).fetchone()
    print(f"  Duplicated seq groups:    {dup_groups[0]:>10,}")
    print(f"  Copies in dup groups:     {dup_groups[1]:>10,}")

    # copy-count distribution
    print(f"\n── Copy-count distribution (for duplicated sequences) ─")
    print(f"  {'copies':>8}  {'# sequences':>12}  {'% of dup groups':>16}")
    rows_dist = con.execute("""
        SELECT copies, COUNT(*) AS n_seqs
        FROM (
            SELECT sequence, COUNT(*) AS copies
            FROM seqs
            GROUP BY sequence
            HAVING COUNT(*) > 1
        )
        GROUP BY copies
        ORDER BY copies
        LIMIT 20
    """).fetchall()
    for copies, n in rows_dist:
        print(f"  {copies:>8}  {n:>12,}  {n/dup_groups[0]*100:>15.1f}%")

    # cross-species duplicates
    print(f"\n── Cross-species duplicates ───────────────────────────")
    cross = con.execute("""
        SELECT COUNT(*) AS n_groups
        FROM (
            SELECT sequence, COUNT(DISTINCT species) AS n_species
            FROM seqs
            GROUP BY sequence
            HAVING COUNT(*) > 1 AND COUNT(DISTINCT species) > 1
        )
    """).fetchone()[0]
    within = dup_groups[0] - cross
    print(f"  Within-species groups:    {within:>10,}")
    print(f"  Cross-species groups:     {cross:>10,}")

    # cross-genus duplicates
    cross_genus = con.execute("""
        SELECT COUNT(*) AS n_groups
        FROM (
            SELECT sequence, COUNT(DISTINCT genus) AS n_genera
            FROM seqs
            GROUP BY sequence
            HAVING COUNT(*) > 1 AND COUNT(DISTINCT genus) > 1
        )
    """).fetchone()[0]
    print(f"  Cross-genus groups:       {cross_genus:>10,}")

    # top 20 most duplicated sequences
    print(f"\n── Top 20 most-copied sequences ───────────────────────")
    print(f"  {'copies':>6}  {'species (distinct)':>25}  {'genera (distinct)':>18}  first species")
    top = con.execute("""
        SELECT
            COUNT(*)                   AS copies,
            COUNT(DISTINCT species)    AS n_species,
            COUNT(DISTINCT genus)      AS n_genera,
            MIN(species)               AS example_species,
            sequence
        FROM seqs
        GROUP BY sequence
        ORDER BY copies DESC
        LIMIT 20
    """).fetchall()
    for copies, n_sp, n_gen, ex_sp, _seq in top:
        print(f"  {copies:>6}  {n_sp:>25}  {n_gen:>18}  {ex_sp}")

    # species with the most redundant sequences
    print(f"\n── Top 20 species by duplicate copies ─────────────────")
    print(f"  {'total':>7}  {'unique':>7}  {'dups':>7}  species")
    top_sp = con.execute("""
        SELECT
            species,
            COUNT(*)               AS total,
            COUNT(DISTINCT sequence) AS unique_seqs,
            COUNT(*) - COUNT(DISTINCT sequence) AS dup_copies
        FROM seqs
        GROUP BY species
        HAVING COUNT(*) > COUNT(DISTINCT sequence)
        ORDER BY dup_copies DESC
        LIMIT 20
    """).fetchall()
    for sp, total_sp, uniq_sp, dup_sp in top_sp:
        print(f"  {total_sp:>7}  {uniq_sp:>7}  {dup_sp:>7}  {sp}")

    if args.write_deduped:
        dedup_path = Path(args.write_deduped)
        dedup_path.parent.mkdir(parents=True, exist_ok=True)
        seen = set()
        written = skipped = 0
        with open(dedup_path, 'w') as out:
            for species, seq, raw_header in raw_records:
                key = (species, seq)
                if key in seen:
                    skipped += 1
                    continue
                seen.add(key)
                out.write(f"{raw_header}\n{seq}\n")
                written += 1
        print(f"\n── Deduped FASTA ──────────────────────────────────────")
        print(f"  Written:  {written:,}  (cross-species duplicates preserved)")
        print(f"  Removed:  {skipped:,}  within-species duplicate copies")
        print(f"  Output:   {dedup_path}")


if __name__ == '__main__':
    main()
