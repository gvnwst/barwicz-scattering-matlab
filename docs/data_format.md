# Data format

Everything the scattering model consumes is in `data/`. There are two kinds of
file, one per pipeline stage boundary.

All files are MATLAB v7 `.mat` (zlib-compressed), readable by MATLAB, GNU
Octave, and `scipy.io.loadmat` alike.

---

## 1. Boundary-field files — `data/boundary_fields/`

Stage-1 output, stage-2 input. **This is the format to match if you want to run
the model on your own waveguides.**

| file | system | n_core | grid | cells | size |
|---|---|---|---|---|---|
| `boundary_fields_SiON_v3_e.mat` | SiON/SiO₂ | 1.50 | 70 × 52 | 2711 | 339 KB |
| `boundary_fields_Si3N4_v3_d.mat` | Si₃N₄/SiO₂ | 2.00 | 40 × 41 | 1130 | 137 KB |
| `boundary_fields_Si_v3.mat` | Si/SiO₂ | 3.50 | 32 × 19 | 353 | 45 KB |

Each holds one variable, `data`, a 2-D struct array indexed `data(ia, ib)`:
`ia` indexes core **height**, `ib` indexes core **width**. The populated region
is triangular, not rectangular — each row runs from a minimum width up to the
single-mode cutoff for that height, so most rows are shorter than the array.
A cell is unpopulated when `data(ia,ib).a` is empty or zero, and that is exactly
the test the downstream loops use to find the end of a row.

### What is in a cell

The striking thing about this interface is how little of the mode survives it.
No field maps are stored — each mode is reduced to **three numbers**, all of them
integrals of the modal field along the two rough sidewalls. A whole 2711-cell
sweep is 339 KB.

```
data(ia,ib)
  .dxy              1x2   grid spacing [um], same for every cell
  .P                      the cross-section that was solved
      .l0                 free-space wavelength [m]          (1.55e-6)
      .nco                core refractive index
      .ncl                cladding refractive index
      .wd                 core width  [m]
      .t                  core height [m]
      .MX, .MY            mode-solver grid size
  .a                      core height [m]   (same as P.t)
  .b                      core width  [m]   (same as P.wd)
  .beta             2x1   propagation constants of both solved modes [rad/m]
  .TE                     the TE-like (x-polarized) mode
      .beta               its propagation constant [rad/m]
      .Phi_squared        sidewall field intensity / guided power  [1/m]
      .R_YZ               |Ez|^2 / |Ex|^2 integrated on the sidewall
      .Ex, .F, .P         zeroed on save (see below)
  .TM                     the TM-like (y-polarized) mode
      .beta               its propagation constant [rad/m]
      .Phi_squared        sidewall field intensity / guided power  [1/m]
      .R_XZ               |Ez|^2 / |Ey|^2 integrated on the sidewall
      .Ey, .F, .P         zeroed on save
  .x_bound          1x2   grid indices of the two sidewalls
  .y_bound          1x2   grid indices of the top and bottom walls
```

`.TE.F`, `.TE.P`, `.TE.Ex` (and the TM equivalents) are scratch variables that
the stage-1 script sets to `0` before saving — they held the full field arrays
during the sweep. They are present in every file as scalar zeros. Do not try to
read a field profile out of them.

### The three quantities that matter

For each polarization, stage 2 uses only `beta`, `Phi_squared` and `R_YZ`/`R_XZ`,
plus the geometry and material indices from `P`. Writing `s` for a sidewall and
integrating vertically over the core height:

- **`Phi_squared`** — total electric-field intensity on one sidewall, normalized
  to the guided power:

  ```
  Phi_squared = ( ∫|Ex|² + ∫|Ey|² + ∫|Ez|² )  /  2P
  ```

  summed over both sidewalls and divided by 2, i.e. the average over the two
  interfaces. This is the strength of the equivalent volume current source.
  `P` is the guided power of that same mode, so the ratio is independent of
  field normalization.

- **`R_YZ`** (TE) and **`R_XZ`** (TM) — the longitudinal-to-transverse intensity
  ratio on the sidewall, `∫|Ez|² / ∫|Ex|²` and `∫|Ez|² / ∫|Ey|²`. These set the
  polarization mix of the radiating current and hence the radiation pattern.

- **`beta`** — the real propagation constant, which fixes which roughness
  spatial frequencies are phase-matched to radiation.

All three come from the sidewall integrals in `boundary_fields_*.m`; the
corresponding lines are the ones computing `TE.R_YZ`, `TM.R_XZ`,
`TE.Phi_squared` and `TM.Phi_squared`.

### Coordinate convention

The mode solver these files came from used a left-handed frame inherited from a
cylindrical code path, and its x/y are transposed relative to the Payne–Lacey
convention used in the scattering literature — the stage-1 script says as much
in a comment. In the stored data:

- **x** is horizontal, across the waveguide width `b`. The rough sidewalls are
  the two constant-x planes. The TE-like mode is the one with dominant `Ex`.
- **y** is vertical, along the waveguide height `a`. Top and bottom walls are
  assumed smooth.
- **z** is the propagation direction.

Because only `|E|²` and real powers survive into the stored quantities, the
handedness cancels; see `modesolve_cross_section.m` for the full discussion.

### Roughness parameters are *not* in this file

`Phi_squared`, `R_*` and `beta` describe the waveguide. The roughness statistics
(correlation length `Lc` and variance `sigma²`) enter in stage 2, where they are
hard-coded as `Lc = 50*nm` and `sigma = 1*nm` near the top of
`scat_hic_integral_final_v3.m`. The published results all assume an exponential
autocorrelation model with `Lc = 50 nm`. Because the leading-order theory is
linear in `sigma²`, the stage-2 output is reported per unit variance and you
rescale afterwards — but changing `Lc` means rerunning stage 2.

---

## 2. Integrated files — `data/integrated/`

Stage-2 output, stage-3 input. Named `<boundary-field name>_integrated.mat`.

| file | size |
|---|---|
| `boundary_fields_SiON_v3_e_integrated.mat` | 778 KB |
| `boundary_fields_Si3N4_v3_d_integrated.mat` | 329 KB |
| `boundary_fields_Si_v3_integrated.mat` | 109 KB |

These are whole-workspace dumps (44 variables) written by
`scat_hic_integral_final_v3.m`, so they carry the original `data` struct array
along with the computed results. Each result is a matrix on the same `(ia, ib)`
grid. The ones that matter:

| variable | meaning |
|---|---|
| `PowerHIC_TE_total`, `PowerHIC_TM_total` | **scattering loss in (dB/cm)/nm²** — the quantity plotted in Fig. 12. Multiply by the expected roughness variance in nm² to get dB/cm. |
| `PowerHIC_TE_factor`, `PowerHIC_TM_factor` | **radiation efficiency** — the 3-D loss divided by the 2-D Payne–Lacey loss, the quantity plotted in Fig. 11. Values below 1 mean the 2-D theory overestimates. |
| `PaynePower_TE`, `PaynePower_TM` | the 2-D Payne–Lacey reference loss used as that denominator |
| `matrix_a`, `matrix_b` | core height and width [m] for each cell |
| `matrix_phi_squared_TE`, `matrix_phi_squared_TM` | the sidewall overlap carried through from stage 1 |
| `PowerHIC_TE_{Y,Z,YZ}_cos_{1,2}`, `PowerHIC_TM_{X,Z}_{cos,sin}_{1,2}` | the individual radiated-power integrals that are summed into the totals |

Unpopulated cells are zero, and the same `data(ia,ib).a > 0` test identifies the
filled region.

---

## 3. Stage-3 output

Stage 3 (`scat_hic_interpolate_*.m`) maps the triangular `(ia, ib)` results onto
a rectangular height/width mesh, pads the edges, and spline-interpolates onto a
mesh four times finer, blanking out-of-range points with `NaN`. It saves its
whole workspace under its own name. The variables to plot are

```
mesh_a_i, mesh_b_i                        height and width [m]
mesh_PowerHIC_TE_total_i, ..._TM_total_i    Fig. 12 panels
mesh_PowerHIC_TE_factor_i, ..._TM_factor_i  Fig. 11 panels
```

`examples/example3_reproduce_fig11_fig12.m` extracts exactly these into a single
tidy `stage3_meshes.mat` with flat names like `Si_te_total`.
