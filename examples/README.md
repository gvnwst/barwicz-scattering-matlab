# Examples

Run these from this directory, in GNU Octave or MATLAB:

```octave
cd examples
example1_loss_lookup
```

Each one adds `../src/waveguide_scattering` to the path itself and writes its
output under `../output/`, which is gitignored. Nothing here modifies anything
in `data/` or `src/`.

| | what it shows | runtime | needs a mode solver |
|---|---|---|---|
| **1** `example1_loss_lookup.m` | Predicts the scattering loss of one waveguide in dB/cm and checks it against the number quoted in the caption of Fig. 12. Start here. | ~1 min | no |
| **2** `example2_single_cell.m` | Runs the scattering model itself on a few cross-sections: prints the dozen numbers that go in, runs the quadrature, prints the loss that comes out, and diffs it against the shipped results. Read this one to understand the interface. | ~30 s per cell, 4 cells by default | no |
| **3** `example3_reproduce_fig11_fig12.m` | Regenerates all twelve panels of Figs. 11 and 12 from the shipped stage-2 results, for all three material systems and both polarizations. | ~8 min | no |
| **4** `example4_validate_stage1.m` | Checks a mode solver you have plugged into `modesolve_cross_section.m` by regenerating part of the Si sweep and diffing it against the reference data. Tells you what to do instead if you have not. | ~2 min per row | **yes** |

## Plotting

`plot_fig11_fig12.py` renders the publication-quality composites — filled,
labelled contours on Tymon Barwicz's own levels — from the mesh file example 3
writes. It is what produced the images in `../figures/`.

```bash
uv run examples/plot_fig11_fig12.py      # from the repository root
```

Example 3 also attempts its own plots in Octave, as a quick look. Those are
unlabelled line contours: Octave's OpenGL toolkits print black frames when no
display is attached, so the example falls back to the gnuplot toolkit, which in
turn renders filled contours badly across the `NaN` region outside the
single-mode band. If you want the real figures, use the Python script.
