#!/bin/bash
cd ~/pi0eta_ampmodel_fit
for d in runs/*/; do
  t=$(basename "$d")
  [ -f "$d/fitpar.npy" ] || continue
  ~/.venv/bin/python3 plots.py "$t" > /dev/null 2>&1 && echo "replotted $t" || echo "FAILED $t"
done
echo "REPLOT DONE"
