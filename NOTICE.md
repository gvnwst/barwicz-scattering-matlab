# Provenance and attribution

## The model

The scattering-loss model and all of its numerical code are the work of
**Tymon Barwicz**, developed at the Research Laboratory of Electronics, MIT, in
2004–2005, with **Hermann A. Haus** (deceased), and published as

> T. Barwicz and H. A. Haus, "Three-Dimensional Analysis of Scattering Losses
> Due to Sidewall Roughness in Microphotonic Waveguides," *Journal of Lightwave
> Technology*, vol. 23, no. 9, pp. 2719–2732, September 2005.
> doi:10.1109/JLT.2005.850816

Tymon Barwicz provided the source files that this repository is built from.

## This packaging

**Gavin West** identified the working subset of the files provided, ported them
to GNU Octave, regenerated every dataset and figure here, and assembled the
repository. The work is described in [docs/validation.md](docs/validation.md).

The code is deliberately *not* refactored. Apart from the four Octave
compatibility fixes listed in that document — each marked `OCTAVE PORT NOTE`
in the file where it appears, and none of which changes a computed value. The `.m` files are as Tymon wrote them.

## Third-party code included

`src/waveguide_scattering/adaptlob.m` — adaptive Lobatto quadrature.
**Walter Gautschi**, bylined 08/03/98 in the file itself, implementing the
algorithm of W. Gander and W. Gautschi, "Adaptive Quadrature — Revisited",
*BIT Numerical Mathematics* 40(1), 84–101 (2000), after Gander,
*Computermathematik*, Birkhäuser, 1992. The paper cites this routine as the
quadrature used for the radiation-profile integrals.

It arrived as part of the original source tree with that attribution line and
no license grant, and is redistributed here on the same basis — **the MIT
license in [LICENSE](LICENSE) does not cover it.**

It is not quite the published version: Tymon added a `global termination2` flag so that the "Interval contains no more machine number" warning fires once per call rather than on every terminal
subinterval, and a `tol < eps` clamp. This doesn't affect the computation. Nothing
else in the file is local work.

## Third-party code not included

The 2-D full-vector finite-difference mode solver used to produce the modal
fields for the paper was written by **Milos Popovic** and is not included here.

The two call sites in the stage-1 sweep drivers were replaced by stubs with
neutral names, `modesolve_cross_section.m` and `mode_power.m`, which document
the required behavior without reproducing the original interface. See
[docs/modesolver_interface.md](docs/modesolver_interface.md).

The three `boundary_fields_*.m` sweep drivers in `src/` were themselves derived
by Tymon from a mode-solver driver script he was given; the physics in
them — the geometry sweep, the TE/TM sorting, the sidewall localization and the
boundary-field reduction — is his. Their original header, which credited the
driver script they started from, has been replaced by a provenance note, and
they do not run as shipped.

## Data

The `.mat` files in `data/` were regenerated from scratch in August 2026 by
Gavin West under GNU Octave 11.3, using the code in `src/` together with the
mode solver that is not distributed here. They are not copies of Tymon
Barwicz's original result files, though they agree with them — see
[docs/validation.md](docs/validation.md) for the comparison.

The figures in `figures/` were rendered from those datasets by
`examples/plot_fig11_fig12.py`. They are reproductions of Figs. 11 and 12 of
the paper, not reproductions *of the published images*; the contour levels are
Tymon Barwicz's, taken from his plotting scripts.

## License

MIT, © Tymon Barwicz and Gavin West — see [LICENSE](LICENSE), which also
records the one carve-out, `adaptlob.m`, described above.
