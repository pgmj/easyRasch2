#!/bin/zsh
# Overnight chain: wait for the mechanism study to finish, then extend the
# decision-stability study from 3 to 6 datasets. Detached; survives the session.
cd "$(dirname "$0")" || exit 1

echo "[$(date)] chain started; waiting for infit_misfit_mechanism.pdf"

# 1. wait for the mechanism render to produce its PDF
while [ ! -f infit_misfit_mechanism.pdf ]; do
  if grep -qiE "Error in" render_mech.log 2>/dev/null; then
    echo "[$(date)] mechanism render errored; continuing anyway"
    break
  fi
  sleep 120
done

# 2. wait for its workers to exit
while [ "$(ps -eo command | grep -c '[r]md.R')" -gt 0 ]; do
  sleep 60
done
echo "[$(date)] mechanism study done; backing up stability cache"
cp infit_decision_stability.rds infit_decision_stability.rds.bak1200 2>/dev/null

echo "[$(date)] launching stability extension (3 -> 6 datasets)"
quarto render infit_decision_stability.qmd > render_stability_6.log 2>&1
echo "[$(date)] stability render exited with status $?"
