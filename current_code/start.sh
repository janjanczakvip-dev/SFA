#!/bin/bash
set -euo pipefail
cd "$(dirname "$0")"

export OMP_NUM_THREADS="${OMP_NUM_THREADS:-32}"

# 1 for 1D, 2 for Hydrogen 3D
model=1
Ipion=1.0          # CLI placeholder; the Volkov phase uses Ip = 0.5
Ncycles=5
omega=0.06         # laser frequency in a.u.
F=0.06             # laser amplitude in a.u.
pmax=4.0
dp=0.001
outfile="spectrum.dat"

echo "Compiling GL-panel..."
make

echo "Running: model=$model nc=$Ncycles omega=$omega F=$F pmax=$pmax dp=$dp threads=$OMP_NUM_THREADS"
./GL-panel "$Ipion" "$Ncycles" "$omega" "$F" "$model" "$pmax" "$dp" "$outfile"

if command -v python3 >/dev/null 2>&1; then
    python3 plot_spectrum.py "$outfile"
fi
