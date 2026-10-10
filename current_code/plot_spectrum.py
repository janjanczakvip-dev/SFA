#!/usr/bin/env python3
"""Plot a GL-panel spectrum file (column 1 = p*0.55, column 2 = |integral|^2)."""

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

if len(sys.argv) < 2:
    print("Usage: python3 plot_spectrum.py spectrum.dat [out_prefix]")
    sys.exit(1)

path = Path(sys.argv[1])
out = Path(sys.argv[2]) if len(sys.argv) > 2 else path.with_suffix("")

data = np.loadtxt(path)
p = data[:, 0] / 0.55
y = data[:, 1]

fig, axes = plt.subplots(2, 1, figsize=(8.5, 7.2), sharex=True)
axes[0].semilogy(p, np.maximum(y, 1e-30), color="C3", lw=1.2)
axes[0].set_ylabel(r"$|\int|^2$")
axes[0].set_title(path.name)
axes[0].grid(True, which="both", alpha=0.3)

axes[1].plot(p, y, color="C3", lw=1.2)
axes[1].set_xlabel(r"$p$")
axes[1].set_ylabel(r"$|\int|^2$")
axes[1].grid(True, alpha=0.3)
fig.tight_layout()
fig.savefig(f"{out}.png", dpi=150)
fig.savefig(f"{out}.pdf")
print(f"Wrote {out}.png and {out}.pdf")
