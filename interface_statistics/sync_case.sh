#!/usr/bin/env bash
# Copy the small statistics products of a simulation case into this repo.
#
#   ./sync_case.sh /scratch/gpfs/DEIKE/<user>/multiphase_cases/<case>  [name]
#
# Only text tables, the per-rank accumulator state, the generated PNGs, and the
# logs are copied. Dumps, restarts, fields, slices and PPMs are never touched.
set -euo pipefail

SRC=${1:?usage: sync_case.sh <case-dir> [name]}
NAME=${2:-$(basename "$SRC")}
DEST="$(cd "$(dirname "$0")" && pwd)/data/$NAME"

[ -d "$SRC/statistics" ] || { echo "no statistics/ in $SRC" >&2; exit 1; }

mkdir -p "$DEST/statistics/plots" "$DEST/logs"

cp "$SRC"/statistics/eta_stats_window.out   "$DEST/statistics/"
cp "$SRC"/statistics/eta_pdf_window_*.out   "$DEST/statistics/"
cp "$SRC"/statistics/eta_stats_local_*.bin  "$DEST/statistics/" 2>/dev/null || true
cp "$SRC"/statistics/eta_stats_state.bin    "$DEST/statistics/" 2>/dev/null || true

# Skip the eta_pdf_latest_* symlinks; git would store them as links to nothing.
for f in "$SRC"/statistics/plots/*.png; do
    [ -e "$f" ] || continue
    [ -L "$f" ] || cp "$f" "$DEST/statistics/plots/"
done

cp "$SRC"/log_*.out "$DEST/logs/" 2>/dev/null || true
cp "$SRC"/run_*.sh  "$DEST/"      2>/dev/null || true

# out.log is huge; keep only the parameter header, which carries T0 = Tp.
if [ -f "$SRC/out.log" ]; then
    sed -n '1,82p' "$SRC/out.log" > "$DEST/case_params.txt"
fi

echo "synced $NAME -> $DEST  ($(du -sh "$DEST" | cut -f1))"
