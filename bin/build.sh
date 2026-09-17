#!/usr/bin/env bash
# Download all reference database sources and convert them to normalised FASTA.
#
# Reads datasets.yaml (repo root) for source URLs, extraction steps, and
# conversion commands. Each dataset is processed independently; pass one or
# more short_names as arguments to run only those datasets.
#
# Usage:
#   bash bin/build.sh                                        # download + convert all
#   bash bin/build.sh gtdb pr2                               # selected datasets
#   bash bin/build.sh --download-only                        # download + prepare only
#   bash bin/build.sh --convert-only                         # convert only (skip download)
#   bash bin/build.sh --skip-udb                             # skip the final UDB build
#   bash bin/build.sh --convert-only gtdb pr2                # flags and filters can combine
#   bash bin/build.sh --list                                 # print available datasets
#   bash bin/build.sh --config small12s.yaml                 # use a custom config file
#   bash bin/build.sh --output-name small_12s                # set output FASTA/UDB base name
#   bash bin/build.sh --source-dir /path/to/storage           # store source data on external storage
#   bash bin/build.sh --output-dir /path/to/storage           # write FASTAs and UDB to external storage
#
# Requirements (cross-platform):
#   yq       mikefarah v4 — NOT the Python jq-wrapper 'yq' (https://github.com/mikefarah/yq)
#   vsearch  only for the final UDB build (https://github.com/torognes/vsearch)
#   python3  >= 3.9, plus per-converter deps (duckdb, pandas, openpyxl)
#   curl, unzip, tar

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
CONFIG="$REPO_ROOT/datasets.yaml"
cd "$REPO_ROOT"

# ── dependency check ──────────────────────────────────────────────────────────
require_cmd() {  # require_cmd <command> <install hint>
    if ! command -v "$1" &>/dev/null; then
        echo "Error: '$1' is required but was not found on PATH." >&2
        echo "  Install: $2" >&2
        exit 1
    fi
}

# ── UDB verification helper (Tier 3) ─────────────────────────────────────────
# Emit up to <n> records from the start and <n> from the end of a FASTA, without
# reading the whole file: awk exits after n records, and tail seeks. Sampling the
# per-dataset parts rather than the combined FASTA guarantees every dataset is
# represented — including ones far smaller than an even step across the combined
# file would ever land on — and, because the parts are concatenated in order, the
# samples still span the whole index positionally.
sample_records() {  # sample_records <fasta> <n>
    local f="$1" n="$2"
    [[ -s "$f" ]] || return 0
    awk -v n="$n" '/^>/ { c++ } c > n { exit } { print }' "$f"
    tail -c 262144 "$f" | awk -v n="$n" '
        /^>/ { started = 1; c++ }
        !started { next }
        c > n { exit }
        { print }'
}

# ── download verification (Tier 0) ────────────────────────────────────────────
# Runs on every downloaded file, and on every cached file before it is trusted.
# Catches the two failure modes that otherwise reach a shipped index silently:
# an endpoint returning an HTML error page, and a truncated archive.
verify_download() {  # verify_download <path>
    local f="$1" kind
    if [[ ! -s "$f" ]]; then
        echo "    verify: FAILED — empty file" >&2
        return 1
    fi
    kind=$(file -b "$f" 2>/dev/null || echo unknown)
    case "$kind" in
        *HTML*|*"XML document"*)
            echo "    verify: FAILED — got '$kind', i.e. an error page rather than data" >&2
            return 1 ;;
    esac
    case "$f" in
        *.zip)
            unzip -qt "$f" >/dev/null 2>&1 || { echo "    verify: FAILED — corrupt or truncated zip" >&2; return 1; } ;;
        *.tgz|*.tar.gz)
            tar -tzf "$f" >/dev/null 2>&1  || { echo "    verify: FAILED — corrupt or truncated tar.gz" >&2; return 1; } ;;
        *.gz)
            gzip -t "$f" 2>/dev/null       || { echo "    verify: FAILED — corrupt or truncated gzip" >&2; return 1; } ;;
    esac
    return 0
}

require_cmd yq      "mikefarah yq v4 — https://github.com/mikefarah/yq  (macOS: brew install yq; Linux: snap install yq, or download the binary)"
require_cmd curl    "your OS package manager  (Debian/Ubuntu: apt install curl; Fedora/RHEL: dnf install curl; macOS: brew install curl)"
require_cmd python3 "https://www.python.org  (Debian/Ubuntu: apt install python3; macOS: brew install python)"

# yq must be the mikefarah Go implementation v4. The Python 'yq' (a jq wrapper)
# and mikefarah v3 use incompatible syntax and would fail cryptically downstream.
yq_version="$(yq --version 2>&1 || true)"
if ! printf '%s\n' "$yq_version" | grep -Eq 'mikefarah|version v?4'; then
    echo "Error: the installed 'yq' is not the required mikefarah v4 implementation." >&2
    echo "  Found:    $yq_version" >&2
    echo "  Required: mikefarah yq v4 — https://github.com/mikefarah/yq" >&2
    echo "  (The Python 'yq' jq-wrapper and mikefarah v3 are NOT compatible.)" >&2
    exit 1
fi

# ── argument parsing ──────────────────────────────────────────────────────────
DO_DOWNLOAD=true
DO_CONVERT=true
DO_UDB=true
REQUESTED=""   # colon-delimited list of requested short_names, empty = all
OUTPUT_NAME="gbif_dna_taxonomy_annotation"
SOURCE_DIR="$REPO_ROOT/source-data"
OUTPUT_DIR="$REPO_ROOT/output/fasta"
DATASET_FASTAS=()  # FASTAs produced by this run, in order
DO_LIST=false
CURL_PROGRESS="--no-progress-meter"
[[ -t 1 ]] && CURL_PROGRESS="--progress-bar"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --download-only) DO_CONVERT=false ;;
        --convert-only)  DO_DOWNLOAD=false ;;
        --skip-udb)      DO_UDB=false ;;
        --list|--help)   DO_LIST=true ;;  # handled after parsing (args are consumed by this loop)
        --config)
            shift
            [[ $# -eq 0 ]] && { echo "Error: --config requires a path argument" >&2; exit 1; }
            CONFIG="$1"
            # Resolve relative paths against repo root
            [[ "$CONFIG" != /* ]] && CONFIG="$REPO_ROOT/$CONFIG"
            [[ -f "$CONFIG" ]] || { echo "Error: config file not found: $CONFIG" >&2; exit 1; }
            ;;
        --output-name)
            shift
            [[ $# -eq 0 ]] && { echo "Error: --output-name requires a name argument" >&2; exit 1; }
            OUTPUT_NAME="$1"
            ;;
        --source-dir)
            shift
            [[ $# -eq 0 ]] && { echo "Error: --source-dir requires a path argument" >&2; exit 1; }
            SOURCE_DIR="$1"
            [[ "$SOURCE_DIR" != /* ]] && SOURCE_DIR="$REPO_ROOT/$SOURCE_DIR"
            ;;
        --output-dir)
            shift
            [[ $# -eq 0 ]] && { echo "Error: --output-dir requires a path argument" >&2; exit 1; }
            OUTPUT_DIR="$1"
            [[ "$OUTPUT_DIR" != /* ]] && OUTPUT_DIR="$REPO_ROOT/$OUTPUT_DIR"
            ;;
        -*) echo "Unknown flag: $1" >&2; exit 1 ;;
        *)  REQUESTED="$REQUESTED:$1:" ;;
    esac
    shift
done

# ── helpers ───────────────────────────────────────────────────────────────────
count=$(yq '.datasets | length' "$CONFIG")

list_datasets() {
    for i in $(seq 0 $((count - 1))); do
        sn=$(yq ".datasets[$i].short_name" "$CONFIG")
        tg=$(yq ".datasets[$i].target_gene" "$CONFIG")
        ver=$(yq ".datasets[$i].version // \"?\"" "$CONFIG")
        sc=$(yq ".datasets[$i].taxonomic_scope" "$CONFIG")
        printf "  %-16s %-6s %-10s %s\n" "$sn" "$tg" "$ver" "$sc"
    done
}

if [[ "$DO_LIST" == true ]]; then
    echo "Available datasets:"
    list_datasets
    exit 0
fi

# vsearch is only needed for the final UDB build — check it before doing any work.
if [[ "$DO_UDB" == true && "$DO_CONVERT" == true ]]; then
    require_cmd vsearch "https://github.com/torognes/vsearch  (conda: conda install -c bioconda vsearch; Debian/Ubuntu: apt install vsearch; macOS: brew install vsearch)"
fi

# ── main loop ─────────────────────────────────────────────────────────────────
for i in $(seq 0 $((count - 1))); do
    short_name=$(yq ".datasets[$i].short_name" "$CONFIG")
    target_gene=$(yq ".datasets[$i].target_gene" "$CONFIG")

    # Skip if a filter was given and this dataset isn't in it
    if [[ -n "$REQUESTED" && "$REQUESTED" != *":$short_name:"* ]]; then
        continue
    fi

    echo ""
    echo "════════════════════════════════════════"
    echo "  $short_name  ($target_gene)"
    echo "════════════════════════════════════════"

    dir="$SOURCE_DIR/$short_name"
    mkdir -p "$dir"

    # ── download endpoints ──────────────────────────────────────────────────
    if [[ "$DO_DOWNLOAD" == true ]]; then
        endpoint_count=$(yq ".datasets[$i].endpoints | length" "$CONFIG")
        curl_flags=$(yq ".datasets[$i].curl_flags // \"\"" "$CONFIG")

        if [[ "$endpoint_count" -eq 0 ]]; then
            echo "  No endpoints configured — skipping download."
        else
            for j in $(seq 0 $((endpoint_count - 1))); do
                url=$(yq ".datasets[$i].endpoints[$j]" "$CONFIG")
                filename=$(basename "$url")
                dest="$dir/$filename"

                if [[ -f "$dest" ]]; then
                    if verify_download "$dest"; then
                        echo "  Already present: $filename"
                    else
                        echo "  Cached $filename failed verification — re-downloading" >&2
                        rm -f "$dest"
                    fi
                fi

                if [[ ! -f "$dest" ]]; then
                    echo "  Downloading $filename …"
                    # Download to .part and rename only after verifying, so an
                    # interrupted transfer is never cached as complete.
                    rm -f "$dest.part"
                    # shellcheck disable=SC2086
                    curl -L --fail --show-error $CURL_PROGRESS \
                         --retry 3 --retry-delay 5 --connect-timeout 30 \
                         -w "    http %{http_code}  %{size_download} bytes  %{time_total}s\n" \
                         $curl_flags -o "$dest.part" "$url"
                    if ! verify_download "$dest.part"; then
                        rm -f "$dest.part"
                        echo "  Aborting: $filename did not verify." >&2
                        exit 1
                    fi
                    mv "$dest.part" "$dest"
                fi
            done
        fi

        # ── prepare (extraction etc.) ─────────────────────────────────────
        prepare_cmd=$(yq ".datasets[$i].prepare_cmd // \"\"" "$CONFIG")
        prepare_cmd="${prepare_cmd//source-data\//$SOURCE_DIR/}"
        prepare_sentinel=$(yq ".datasets[$i].prepare_sentinel // \"\"" "$CONFIG")
        prepare_sentinel="${prepare_sentinel//source-data\//$SOURCE_DIR/}"

        if [[ -n "$prepare_cmd" ]]; then
            if [[ -n "$prepare_sentinel" && -e "$prepare_sentinel" ]]; then
                echo "  Preparation already done (sentinel exists: $prepare_sentinel)"
            else
                echo "  Preparing …"
                eval "$prepare_cmd"
            fi
        fi
    fi

    # ── convert ─────────────────────────────────────────────────────────────
    if [[ "$DO_CONVERT" == true ]]; then
        convert_cmd=$(yq ".datasets[$i].convert_cmd" "$CONFIG")
        convert_cmd="${convert_cmd//source-data\//$SOURCE_DIR/}"
        convert_cmd="$convert_cmd --output-dir \"$OUTPUT_DIR\""
        echo "  Converting …"
        eval "$convert_cmd"
        # Derive the output FASTA path: last positional argument before any -- flags
        fasta_stem=$(echo "$convert_cmd" | sed 's/ --[a-z].*//' | awk '{print $NF}')
        DATASET_FASTAS+=("$OUTPUT_DIR/${fasta_stem}.fasta")

        # ── post-process (e.g. within-species dedup) ──────────────────────
        postprocess_cmd=$(yq ".datasets[$i].postprocess_cmd // \"\"" "$CONFIG")
        postprocess_fasta=$(yq ".datasets[$i].postprocess_fasta // \"\"" "$CONFIG")
        if [[ -n "$postprocess_cmd" ]]; then
            postprocess_cmd="${postprocess_cmd//source-data\//$SOURCE_DIR/}"
            postprocess_cmd="${postprocess_cmd//output\/fasta\//$OUTPUT_DIR/}"
            echo "  Post-processing …"
            eval "$postprocess_cmd"
            if [[ -n "$postprocess_fasta" ]]; then
                # Use the post-processed FASTA in the combined output instead of the raw one.
                # (Index explicitly rather than [-1]; negative subscripts need bash 4.3+, not macOS 3.2.)
                DATASET_FASTAS[$((${#DATASET_FASTAS[@]} - 1))]="$OUTPUT_DIR/${postprocess_fasta}.fasta"
            fi
        fi
    fi

done

# ── concatenate all dataset FASTAs into one combined file ─────────────────────
if [[ "$DO_CONVERT" == true ]]; then
    mkdir -p "$OUTPUT_DIR"
    COMBINED="$OUTPUT_DIR/${OUTPUT_NAME}.fasta"
    echo ""
    echo "════════════════════════════════════════"
    echo "  Concatenating all FASTAs"
    echo "════════════════════════════════════════"

    parts=("${DATASET_FASTAS[@]+"${DATASET_FASTAS[@]}"}")

    if [[ ${#parts[@]} -eq 0 ]]; then
        echo "  No FASTAs were produced — nothing to concatenate." >&2
        exit 1
    fi

    cat "${parts[@]}" > "$COMBINED"
    COUNT=$(grep -c "^>" "$COMBINED")
    echo "  Done — $COUNT sequences written to $COMBINED"

    # ── build vsearch UDB index ───────────────────────────────────────────────
    if [[ "$DO_UDB" == true ]]; then
        echo ""
        echo "════════════════════════════════════════"
        echo "  Building vsearch UDB"
        echo "════════════════════════════════════════"
        UDB="$OUTPUT_DIR/${OUTPUT_NAME}.udb"
        LOG="$OUTPUT_DIR/${OUTPUT_NAME}.log"
        vsearch --makeudb_usearch "$COMBINED" --output "$UDB" --log "$LOG"
        echo "  Done — $UDB"

        # ── self-hit smoke test (Tier 3) ─────────────────────────────────────
        # Every count-based check above can pass on a corrupt UDB. This one
        # cannot: sequences taken from the combined FASTA must find themselves
        # in the index that was just built from it.
        echo ""
        echo "  Verifying UDB — self-hit query …"
        SAMPLE="$OUTPUT_DIR/.udb_selftest.fasta"
        HITS="$OUTPUT_DIR/.udb_selftest.tsv"
        PER_PART=4
        : > "$SAMPLE"
        for part in "${parts[@]}"; do
            sample_records "$part" "$PER_PART" >> "$SAMPLE"
        done
        # Key on the WHOLE header, which is what vsearch reports in blast6out
        # column 1 (headers contain no whitespace, so nothing is truncated).
        # Keying on the first field instead would merge distinct records: one
        # GenBank accession appears across several MIDORI gene sets, since a
        # mitochondrial genome contributes a region per gene. sort -u also
        # absorbs a genuine duplicate if a part is smaller than the tail window.
        got=$(awk '/^>/ { print substr($0, 2) }' "$SAMPLE" | sort -u | wc -l | tr -d ' ')
        echo "  Sampled $got sequences from ${#parts[@]} datasets"

        vsearch --usearch_global "$SAMPLE" --db "$UDB" \
                --id 0.99 --maxaccepts 1 --maxhits 1 --maxrejects 32 \
                --blast6out "$HITS" --quiet --threads 4

        matched=$(awk -F'\t' '$3 >= 99.0 {print $1}' "$HITS" | sort -u | wc -l | tr -d ' ')
        if [[ "$got" -eq 0 ]]; then
            echo "  UDB self-hit test FAILED — could not sample the combined FASTA." >&2
            exit 1
        fi
        pct=$(( matched * 100 / got ))
        echo "  Self-hit: $matched/$got ($pct%) sampled sequences found themselves at >=99% identity"
        if [[ "$pct" -lt 95 ]]; then
            echo "  UDB self-hit test FAILED (<95%) — the index is incomplete or corrupt." >&2
            echo "  Leaving $SAMPLE and $HITS in place for inspection." >&2
            exit 1
        fi
        rm -f "$SAMPLE" "$HITS"
        echo "  UDB verified."
    fi
fi

echo ""
echo "Done."
