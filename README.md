# sequence-annotation-db-build

Builds a combined [vsearch](https://github.com/torognes/vsearch) UDB reference database from
public DNA reference sources.

The build downloads each source listed in `datasets.yaml`, converts it to a common FASTA header
format, concatenates everything, and indexes the result. One command, one config file, two
artefacts: a combined FASTA and the UDB built from it.

---

## Requirements

| Tool | Notes |
|---|---|
| `yq` | **mikefarah v4** — *not* the Python jq-wrapper of the same name. The build checks this and fails fast. `brew install yq` / `snap install yq` |
| `vsearch` | Only needed for the final UDB build. `conda install -c bioconda vsearch`, `apt install vsearch`, `brew install vsearch` |
| `curl`, `unzip`, `tar` | Usually present |
| `python3` | >= 3.9 |
| Python packages | `duckdb`, `pandas`, `openpyxl` — `pip install duckdb pandas openpyxl` |

Linux and macOS both work. No network access is needed beyond the source endpoints and, for
one dataset, the NCBI Entrez API.

## Resource requirements

Measured on a full 26-source build (2026-09):

| | |
|---|---|
| Wall time | ~60 min |
| Downloaded sources | ~11 GB |
| Combined FASTA | ~4.7 GB |
| UDB | ~18 GB |
| Peak RAM | ~8 GB (the vsearch index build) |
| **Free disk needed** | **~35 GB**, sources and outputs combined |

Sources and outputs can live on separate volumes — see `--source-dir` and `--output-dir`.

## Quick start

```bash
git clone <this repo> && cd sequence-annotation-db-build

cp .env.example .env          # optional, see Secrets below
bash bin/build.sh --list      # show the configured sources
bash bin/build.sh             # full build
```

Outputs land in `output/fasta/`:

```
gbif_dna_taxonomy_annotation.fasta   combined, all sources
gbif_dna_taxonomy_annotation.udb     the vsearch index
gbif_dna_taxonomy_annotation.log     vsearch's own build log
<source>.fasta                       one per source, kept for inspection
```

Both the name and the location are configurable (`--output-name`, `--output-dir`).

## Options

```
bash bin/build.sh                              full build
bash bin/build.sh gtdb pr2                     only the named sources  (see caveat below)
bash bin/build.sh --list                       list configured sources and exit
bash bin/build.sh --download-only              fetch and extract, no conversion
bash bin/build.sh --convert-only               convert from already-downloaded sources
bash bin/build.sh --skip-udb                   build the combined FASTA but not the index
bash bin/build.sh --config other.yaml          use a different config file
bash bin/build.sh --output-name small_12s      name the combined FASTA/UDB
bash bin/build.sh --source-dir  /mnt/data      where downloads are stored
bash bin/build.sh --output-dir  /mnt/data/out  where FASTAs and the UDB are written
```

Flags and source filters combine.

Downloads are cached: a source already present and passing verification is not re-fetched, so a
re-run after a failure resumes rather than starting over.

## Secrets

Only one source needs credentials. `ncbi_matk` queries the NCBI Entrez API live at convert time
and has no static download.

```bash
cp .env.example .env     # then fill in NCBI_API_KEY
set -a; source .env; set +a
```

The build runs without a key at NCBI's anonymous rate limit of 3 requests/sec; with one it gets
10/sec. A key is free from https://www.ncbi.nlm.nih.gov/account/ → Settings → API Key Management.

`.env` is gitignored. Never commit it.

## Configuration

`datasets.yaml` is the whole configuration. Each entry describes one source:

```yaml
- short_name: gtdb                      # identifier; names the output FASTA and source dir
  version: R220                         # informational only, not verified — see Caveats
  target_gene: SSU_rRNA_16S_prokaryotic # written verbatim into every header of this source
  taxonomic_scope: Bacteria and Archaea # informational
  citation: >
    Parks DH et al. (2022) ...
  endpoints:                            # downloaded to <source-dir>/<short_name>/
    - https://data.gtdb.ecogenomic.org/...
  curl_flags: --insecure                # optional, appended to the curl call
  prepare_cmd: tar -xzf ...             # optional, run after download (extraction, renames)
  prepare_sentinel: source-data/...     # optional, skip prepare_cmd if this path exists
  convert_cmd: python3 converters/gtdb_to_fasta.py ...   # required
  postprocess_cmd: python3 converters/duplicate_analysis.py ...  # optional
  postprocess_fasta: gtdb_16s_deduped   # optional, use this FASTA in the combined output
```

In `prepare_cmd`, `convert_cmd` and `postprocess_cmd`, the literal prefixes `source-data/` and
`output/fasta/` are rewritten to the active `--source-dir` and `--output-dir`. Write them as
those literals and the paths follow wherever the build is pointed.

`target_gene` is copied into the header of every sequence from that source and is what
downstream consumers filter on, so keep the values consistent across sources.

### Adding a source

Add an entry, and point `convert_cmd` at a script in `converters/` that emits the shared header
format. Existing converters cover Darwin Core archives, plain and gzipped FASTA with various
header conventions, Excel, and a live NCBI query — one of them is usually close to what a new
source needs.

### Header format

Every converter emits the same 23 pipe-separated fields:

```
>ID|accessionNumber|scientificName|decimalLatitude|decimalLongitude|typeStatus|catalogueNumber|
 identifiedBy|taxonRank|country|locality|basisOfRecord|higherClassification|dataset|targetGene|
 domain|kingdom|phylum|class|order|family|genus|species
```

Empty fields are allowed; the count is not. Field values must not contain `|` or `>` — the
converters strip both.

## Verification

Two checks run automatically and fail the build:

- **Every download** is checked for being empty, being an HTML error page rather than data, and
  for archive integrity (`unzip -t` / `tar -tzf` / `gzip -t`). Cached files are re-checked
  before reuse and re-fetched if they fail. Downloads land in `.part` files and are renamed only
  after passing, so an interrupted transfer is never cached as complete.
- **The finished UDB** is queried with ~200 sequences sampled from every source. Each must find
  itself at >=99% identity; below 95% the build fails and leaves the sample and results in
  `output/fasta/.udb_selftest.*` for inspection. This is the only check that catches a
  structurally plausible but corrupt index.

## Running on a schedule

- Allow well over an hour. `ncbi_matk` is a live API query whose duration depends on NCBI and
  which produces no progress output while running.
- Log to a file; the progress bar is automatically suppressed when stdout is not a terminal.
- The build is resumable — cached, verified downloads are skipped on a re-run.
- Exit code is non-zero on failure. There is no partial-success mode: see Caveats.

## Caveats

These are known and deliberate to leave visible rather than hide.

- **A filtered run rebuilds the combined FASTA and UDB from only the selected sources.**
  `bash bin/build.sh gtdb` will replace a full combined index with a gtdb-only one. Use
  `--skip-udb`, or a separate `--output-dir`, when building a subset.
- **One failing source aborts the whole run.** There is no per-source error isolation and no
  summary of what succeeded.
- **Outputs are written in place.** A crash or a full disk during concatenation or indexing
  leaves a truncated artefact where a working one was. There is no atomic promote step and no
  free-space preflight.
- **Most endpoint URLs pin a version** (`GB269`, release dates, `v5.1.0.0`). When an upstream
  rotates, the download fails and the entry needs updating. `version:` in the config is a label,
  not something the build verifies. The exception is `boldistilled`, whose endpoint always
  serves the current release; its converter resolves the release from the extracted filenames
  and records it in `<output>.fasta.version`.
- **`curl_flags: --insecure` is set on the MIDORI sources**, disabling TLS verification for
  those downloads.
- **UNITE is fetched from a mirror** because the official URL sits behind a click-through
  agreement — see the comment on that entry.
