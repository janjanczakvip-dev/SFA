# Current SFA2 wavepacket integrator

Standalone copy of the production code. It evaluates the same 2-D integral as
the original `GL-test.cpp`, with analytic \(S_1,S_2,k\), composite
Gauss–Legendre panels, \(v=s^2\) on the first \(v\)-panel, and one pass over
quadrature nodes for the whole \(p\)-grid.

Sources: `GL-panel.cpp`, `common.hpp`.

## Build

Needs `g++` with OpenMP.

```bash
make
export OMP_NUM_THREADS=32    # optional; default is all cores
```

## Run

```bash
./GL-panel  Ip_ion  nc  omega  F0  model  pmax  dp  [outfile]
            [n_u  n_s  n_s_panels  osc_per_panel  use_s2]
```

Example (5-cycle pulse, \(F_0=\omega=0.06\), model 1, \(p=0\ldots 4\)):

```bash
./GL-panel 1.0 5 0.06 0.06 1 4.0 0.001 spectrum.dat
```

Smoke test (9 momenta, should finish in well under a second):

```bash
./GL-panel 1.0 5 0.06 0.06 1 4.0 0.5 test.dat
```

### Arguments

| argument | meaning | typical |
|---|---|---|
| `Ip_ion` | CLI placeholder (the phase uses \(I_p=0.5\)) | `1.0` |
| `nc` | optical cycles in the \(\sin^2\) pulse | `5` |
| `omega` | laser angular frequency (a.u.) | `0.06` |
| `F0` | peak field (a.u.) | `0.06` |
| `model` | `1` = 1-D, `2` = 3-D H-like | `1` |
| `pmax` | last momentum | `4.0` or `12.0` |
| `dp` | momentum step | `0.001` |
| `outfile` | output path (optional) | `spectrum.dat` |

Optional quadrature overrides (defaults are used if omitted):

| argument | default | role |
|---|---|---|
| `n_u` | 8 | GL order per \(u\)-panel |
| `n_s` | 12 | GL order per \(v\)-panel |
| `n_s_panels` | \( \max(24,\mathrm{round}(40 F_0/\omega)) \) | number of \(v\)-panels |
| `osc_per_panel` | 2 | Volkov oscillations per \(u\)-panel at \(p_{\max}\) |
| `use_s2` | 1 | first \(v\)-panel mapped with \(v=s^2\) |

If you omit `n_s_panels`, it follows the field:
\(n_{s,\mathrm{panels}}=\max(24, \mathrm{round}(40 F_0/\omega))\).
Passing it on the command line turns that off.

\(u\)-panel length already shrinks as \(p_{\max}^2/2+I_p\). Do not also
scale it with \(F_0\).

## Output

Two columns, 16-digit:

```
p * 0.55     |integral|^2
```

The factor \(0.55\) is the original plotting convention. To plot versus
\(p\), divide column 1 by 0.55:

```bash
python3 plot_spectrum.py spectrum.dat
```

Stdout includes `TIMING integrate_sec=…` and `QUAD …` (node count, chosen
`n_s_panels`, \(A_0=F_0/\omega\)).

## Notes

- Pulse length \(t_{\mathrm{length}}=n_c\cdot 2\pi/\omega\), integral to
  \(t_{\max}=1.5\,t_{\mathrm{length}}\).
- No Python Gauss–Legendre file is needed; nodes are generated in C++.
- Compile with `-O3`. Without it the same algorithm is several times slower.
