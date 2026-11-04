#!/usr/bin/env bash
## Run every figure script in order. Test mode (simulated data, ~5 min):
##   LOCUSPICK_PAPER_TEST=1 bash run_all.sh
## Real data: run the 00_* steps in README.md first, then  bash run_all.sh
set -euo pipefail
cd "$(dirname "$0")"
mode=$([ "${LOCUSPICK_PAPER_TEST:-0}" = "1" ] && echo results_test || echo results)
mkdir -p "$mode/logs"
for s in 01_fig2_audit.R 02_fig3_shortcuts.R 03_fig4_global_shift.R 04_fig5_robustness.R \
         05_fig6_singlecell.R 06_fig7_case_study.R; do
  echo "== $s"
  if Rscript "$s" > "$mode/logs/${s%.R}.log" 2>&1; then echo "   ok"
  else echo "   FAILED: see $mode/logs/${s%.R}.log"; tail -5 "$mode/logs/${s%.R}.log"; exit 1; fi
done
Rscript -e 'sessionInfo()' > "$mode/logs/sessionInfo.txt"
echo "All figures in $mode/"
