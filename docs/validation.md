# Validation

The code in `src/` is Tymon Barwicz's 2004–2005 MATLAB, run essentially
unmodified. The question this page answers is whether it still computes what it
computed then, and whether what it computes matches the paper.

Work done August 2026 by Gavin West, under GNU Octave 11.3 on Windows.

## Summary

| check | result |
|---|---|
| Stage 2 under Octave vs. Barwicz's own MATLAB-2005 results | ~1e-16 relative |
| Stage 3 under Octave vs. his saved workspace | bit-exact (coarse), ~5e-16 (interpolated) |
| Stage 1 regenerated vs. his mode data | median ~1e-11 |
| Stage 2 on regenerated stage-1 data | median ~1e-11, max ~4e-6 |
| Final maps, fully regenerated vs. his data | worst case ~2.5%, at the degenerate diagonal |
| Figs. 11 and 12, all six panels, all three systems | reproduced |
| Fig. 12 caption worked example (Si, 200×400 nm, TE, σ²=10 nm²) | 6.67 dB/cm vs. "around 7 dB/cm" |

That last row is what `examples/example1_loss_lookup.m` prints.

## What was checked, in order

**1. Does the 2005 code still run?** A single stage-2 integration cell
`(ia,ib) = (1,1)` was recomputed under Octave from Barwicz's own stage-1 mode
data and compared against his saved integrated results: all 14 computed
quantities matched to ~1e-16 relative. That establishes MATLAB-2005 ↔ Octave-11
numerical equivalence for the whole stage-2 chain — the adaptive Lobatto
quadrature, the six Mathematica-exported integrands, and the Payne–Lacey
reference.

**2. Does stage 3 still run?** `scat_hic_interpolate_Si.m` was rerun on his
integrated data and compared against his saved stage-3 workspace: bit-exact on
the coarse meshes, ~5e-16 after spline interpolation.

**3. Can stage 1 be regenerated?** The full 353-cell Si sweep was rerun from
scratch. It matches his mode data to a median ~1e-11. Only four cells differ by
more than 1%, all of them the near-square cross-sections where the two modes are
degenerate — three of which he hand-patched himself.

**4. Does the whole pipeline close?** Stage 2 on the regenerated stage-1 data
matches his integrated results to a median ~1e-11 (max ~4e-6, degenerate cells
excluded), and the final interpolated maps agree to a worst case ~2.5%,
localized at that same degenerate diagonal.

**5. All three material systems.** Stage 1 → 2 → 3 was run end to end for
Si₃N₄/SiO₂ and SiON/SiO₂ as well, reproducing all six panels of Fig. 11 and all
six of Fig. 12. The composites are in [`figures/`](../figures).

## The degenerate cells

Where the core is square, the TE-like and TM-like modes are degenerate and any
superposition of them is also a mode. Which one a solver returns first is
arbitrary and sensitive to rounding, so the TE/TM assignment there is not
reproducible between solvers or even between runs. `boundary_fields_*.m` prints

```
warning: Both modes are either TE or TM. Degeneracy may create wrong Ez pattern.
```

for those cells and falls back to taking the first mode as TE. This is the one
place where regenerated data legitimately differs from the shipped reference,
and it is also the one visible blemish in the reproduced maps — a thin ridge
along the `a = b` diagonal. The paper's own figures show it too.

## Porting notes

Four changes were needed to run the 2005 MATLAB under Octave. All are marked
`OCTAVE PORT NOTE` in the files, and none changes any computed value:

1. `scat_hic_interpolate_Si.m` — removed a stray subfunction `end` (syntax only).
2. `boundary_fields_*.m` — added `clear TE TM` before the mode-sorting branch.
   The previous iteration leaves `TE.F = 0`, and Octave and modern MATLAB both
   refuse the subsequent `TE.F.Ex = ...` struct assignment over a scalar, which
   MATLAB 2005 allowed. `TE` and `TM` are fully reassigned every iteration, so
   behavior is unchanged. (Gavin West applied the same fix in 2019.)
3. `boundary_fields_Si3N4_v3_d.m` — commented out a `load` of a checkpoint file
   that no longer exists. The `_v3_b`/`_c`/`_d` chain is one sweep resumed after
   crashes; the drivers were verified identical by diff.
4. Run wrappers set
   `warning('off','Octave:possible-matlab-short-circuit-operator')`.
   `adaptlob.m` uses `|` and `&` on scalars, which Octave warns about copiously
   but executes exactly as MATLAB does. `adaptlob.m` itself is untouched.

The SiON stage-1 output is named `boundary_fields_SiON_v3_e` rather than
`_v3`, so that it lines up with the filename `scat_hic_interpolate_SiON.m`
expects. The `_v3` through `_v3_e` drivers are one identical sweep, verified by
diff.

## Reproducing the validation

`examples/example2_single_cell.m` is checks 1 and 4 in miniature: it recomputes
stage-2 cells from the shipped stage-1 data and diffs them against the shipped
stage-2 data. It should report a worst disagreement of exactly zero, since the
quadrature is deterministic and both sides run the same code.

`examples/example4_validate_stage1.m` is check 3, and needs a mode solver.
