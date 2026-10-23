# Plugging in a mode solver

## Why there isn't one here

The model needs a 2-D full-vector waveguide mode solver to turn a waveguide
cross-section into the sidewall quantities it integrates. The solver Tymon
Barwicz used for the paper was written by Milos Popovic and is not included
here, and neither is any part of its API.

Two files stand where it used to be called:

| file | replaces | returns |
|---|---|---|
| `src/waveguide_scattering/modesolve_cross_section.m` | the mode solve | `beta`, `F` |
| `src/waveguide_scattering/mode_power.m` | the power normalization | `P`, `SzTM`, `SzTE` |

Both currently raise `waveguide_scattering:noModeSolver`. Implement them to run
stage 1. Their headers carry the full contract; this page is the context around
it.

## You probably don't need one

Stage 1 is the *only* part of the pipeline that touches a mode solver, and its
output for all three material systems of the paper is shipped in
`data/boundary_fields/`. Stages 2 and 3 — the scattering model itself and the
map generation — run from that data with no solver at all, which is what
examples 1 through 3 do.

You need a solver only to run the model on a waveguide geometry, wavelength or
material system that is not in the shipped sweeps.

## What the solver has to produce

The quantities that leave stage 1 are listed in [data_format.md](data_format.md).
Working back from them, a solver must give you, for each of the two
lowest-order modes of the cross-section:

1. the propagation constant `beta`;
2. `Ex`, `Ey`, `Ez` over the cross-section, from which the sidewall intensity
   integrals `Phi_squared`, `R_YZ` and `R_XZ` are formed;
3. the guided power of the same mode, which normalizes those integrals.

That is all. The absolute field normalization is irrelevant — every stored
quantity is a ratio — so you do not need to match any particular convention.

## The two things that are easy to get wrong

**Grid staggering.** The original solver returned the six field components on
mutually offset Yee-type sub-grids of *different sizes*: with a field grid
`MX x MY`, `Ex` and `Hy` are `(MX-1) x (MY-2)`, while `Ey` and `Hx` are
`(MX-2) x (MY-1)`. `mode_power` relies on `Ex`/`Hy` being co-located and
`Ey`/`Hx` being co-located, and the stage-1 sidewall integrals index `Ex`, `Ey`
and `Ez` with the *same* index ranges despite the half-pixel offsets between
them. That last part is an approximation in the original code, not an artifact
of these stubs; it is second-order in the grid spacing and the published numbers
were produced with it. A replacement solver that pre-interpolates everything
onto one common grid will give slightly different stage-1 numbers — usually
better ones, but no longer bit-comparable with `data/boundary_fields/`.

**Power convention.** `mode_power` does not compute `(1/2)Re{E × H*}`. It uses
the unconjugated product with no factor of ½, and a sign flip on the `Ex·Hy`
term that compensates the left-handed frame the original solver worked in. For
lossless bound modes these agree up to an overall real constant, and since `P`
only ever appears as a denominator the constant cancels — but if you mix
conventions between the power and the field integrals, `Phi_squared` picks up a
spurious factor and every predicted loss is wrong by it. Compute `P` from the
same arrays the sidewall integrals use.

## Checking your implementation

`examples/example4_validate_stage1.m` regenerates the first few rows of the Si
sweep with your solver and compares every stored quantity against
`data/boundary_fields/boundary_fields_Si_v3.mat`, which was produced with the
original solver.

A faithful implementation agrees to around 1e-11. Expect genuine disagreement
only at the near-square cross-sections (`a ≈ b`), where the two modes are
degenerate and the TE/TM assignment is arbitrary; the stage-1 script prints a
warning for those cells, and Tymon hand-patched several of them. See
[validation.md](validation.md).

## Running a new geometry

The three `boundary_fields_*.m` scripts are sweep drivers with the geometry
hard-coded near the top — `min_a`/`max_a`, `min_b`/`max_b`, the step `da`/`db`,
the domain size, the grid spacing `dxy`, and the `P.nco`/`P.ncl`/`P.l0`
material block. The single-mode cutoff that ends each row is the fitted
expression in the `b = ...` line. Copy one, edit those, and you have a sweep for
a different material system.

Mind the cost before you do: see the timings in the main
[README](../README.md#cost).
