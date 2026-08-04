#!/bin/zsh
# Wait for the poly_weak 120-rep run to finish, then launch the decision
# stability study. Detached: survives independently of the session.
cd "$(dirname "$0")" || exit 1

TARGET=4200
CACHE=infit_width_iterations_poly_weak.rds

echo "[$(date)] chain started; waiting for $CACHE to reach $TARGET cells"

# 1. wait for all cells to be computed
while true; do
  N=$(Rscript -e 'cat(tryCatch(nrow(readRDS("infit_width_iterations_poly_weak.rds")$grid), error=function(e) 0))' 2>/dev/null)
  case "$N" in (''|*[!0-9]*) N=0 ;; esac
  if [ "$N" -ge "$TARGET" ]; then
    echo "[$(date)] cache reached $N cells"
    break
  fi
  sleep 300
done

# 2. wait for the render process to exit (PDF regeneration, plots, etc.)
while [ "$(ps -eo command | grep -c '[r]md.R')" -gt 0 ]; do
  sleep 60
done
echo "[$(date)] poly_weak render finished; launching stability study"

# 3. launch
quarto render infit_decision_stability.qmd > render_stability.log 2>&1
echo "[$(date)] stability render exited with status $?"
