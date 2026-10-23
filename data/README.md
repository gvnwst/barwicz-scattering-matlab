# Data

Everything needed to reproduce Figs. 11 and 12 of the paper without a mode
solver. 1.7 MB in total.

All files are MATLAB v7 `.mat`, readable by MATLAB, GNU Octave and
`scipy.io.loadmat`. The format is specified in
[../docs/data_format.md](../docs/data_format.md); what follows is just an index.

## `boundary_fields/` — stage-1 output, 521 KB

The waveguide side of the problem: each cross-section reduced to the handful of
sidewall quantities the scattering model integrates. **This is the format to
match if you want to run the model on your own waveguides.**

| file | system | n_core | heights (nm) | widths (nm) | cross-sections |
|---|---|---|---|---|---|
| `boundary_fields_SiON_v3_e.mat` | SiON/SiO₂ | 1.50 | 550–4000 | 500–4000 | 2711 |
| `boundary_fields_Si3N4_v3_d.mat` | Si₃N₄/SiO₂ | 2.00 | 220–1000 | 200–1000 | 1130 |
| `boundary_fields_Si_v3.mat` | Si/SiO₂ | 3.50 | 135–600 | 120–585 | 353 |

All at λ = 1550 nm, cladding n = 1.45, with both the TE-like and TM-like
fundamental modes. Each sweep covers a triangular region: from a minimum width
up to the single-mode cutoff for that height.

Generated in August 2026 under GNU Octave 11.3 with the code in `../src/` plus
a mode solver that is not distributed here — roughly 9 hours of computation.

## `integrated/` — stage-2 output, 1.2 MB

The model's answers for those same cross-sections: scattering loss per unit
roughness variance, and the radiation efficiency relative to a 2-D
Payne–Lacey analysis. Assumes exponentially correlated roughness with
**Lc = 50 nm**, the value used throughout the paper.

| file | from |
|---|---|
| `boundary_fields_SiON_v3_e_integrated.mat` | the SiON sweep above |
| `boundary_fields_Si3N4_v3_d_integrated.mat` | the Si₃N₄ sweep above |
| `boundary_fields_Si_v3_integrated.mat` | the Si sweep above |

These took about six days of CPU time to produce (see the cost table in the
main [README](../README.md#cost)), which is why they are shipped: with them,
every figure in the paper regenerates in minutes.

## Filenames

The `_v3`, `_v3_d`, `_v3_e` suffixes are Tymon Barwicz's own, carried over
because the stage-3 scripts load these names literally. They mark sweeps
resumed after crashes rather than different physics; the drivers in each chain
were verified identical by diff.
