function [beta, F] = modesolve_cross_section(N, k0, mu_guess, OPTS, nmodes, PML_WIDTH, PML_STRENGTH)
% MODESOLVE_CROSS_SECTION  Stage-1 mode-solver plug-in point (NOT IMPLEMENTED).
%
% The stage-1 sweeps (boundary_fields_*.m) were written against a 2D
% full-vector finite-difference mode solver that is not included with this
% project.  Its call site has been replaced by this stub.  Supply your own
% solver here to run stage 1.
%
% You do NOT need a mode solver to use this repository: precomputed stage-1
% output for all three material systems of the paper ships in
% data/boundary_fields/, and stages 2 and 3 consume that directly.
%
% ---------------------------------------------------------------------------
% INPUTS
% ---------------------------------------------------------------------------
%   N             index cross-section, on a grid twice as fine as the field
%                 grid in each direction:
%                   N.n  (2*MX-1) x (2*MY-1) real refractive index
%                   N.x  (2*MX-1) x 1  x coordinates [m]
%                   N.y  (2*MY-1) x 1  y coordinates [m]
%                 The callers build N with index averaging on the half-grid
%                 points; see the "Dielectric constant averaging" block in
%                 boundary_fields_Si_v3.m.  MX and MY are recorded as P.MX,
%                 P.MY and are 341/321 (Si), 356/346 (Si3N4), 361/341 (SiON).
%   k0            free-space wavenumber 2*pi/lambda0 [rad/m]
%   mu_guess      initial guess for the propagation constant [rad/m]
%                 (the callers pass n_core*k0)
%   OPTS          solver options as set by the callers:
%                   OPTS.fieldmode = 'V'  full-vector
%                   OPTS.coordmode = 'C'  Cartesian
%                   OPTS.eigmode   = 'b'  solve for beta at fixed frequency
%   nmodes        number of modes to return (the callers ask for 2)
%   PML_WIDTH     [left right bottom top] absorber widths [m].  The callers
%                 pass [0 0 0 0], i.e. no absorber: the structures of the
%                 paper are buried and the domain is padded instead.
%   PML_STRENGTH  absorber strength; unused when PML_WIDTH is zero.
%
% ---------------------------------------------------------------------------
% OUTPUTS
% ---------------------------------------------------------------------------
%   beta          nmodes x 1 complex propagation constants [rad/m], sorted
%                 lowest-loss first then highest effective index first.  Only
%                 real(beta) feeds the model; imag(beta) is used for a
%                 diagnostic loss printout.
%
%   F             struct of modal fields on a staggered (Yee-type) grid.  The
%                 six components do NOT share one array size -- each sits on
%                 its own sub-grid.  With the field grid MX x MY:
%
%                   F.Ex   (MX-1) x (MY-2) x nmodes     on (F.Rr, F.Zr)
%                   F.Ey   (MX-2) x (MY-1) x nmodes     on (F.Rz, F.Zz)
%                   F.Ez   (MX-2) x (MY-2) x nmodes     on (F.Rz, F.Zr)
%                   F.Hx   (MX-2) x (MY-1) x nmodes     on (F.Rz, F.Zz)
%                   F.Hy   (MX-1) x (MY-2) x nmodes     on (F.Rr, F.Zr)
%                   F.Hz   (MX-1) x (MY-1) x nmodes
%
%                   F.Rr   (MX-1) x 1   x coords of the Ex/Hy sub-grid [m]
%                   F.Zr   (MY-2) x 1   y coords of the Ex/Hy sub-grid [m]
%                   F.Rz   (MX-2) x 1   x coords of the Ey/Hx sub-grid [m]
%                   F.Zz   (MY-1) x 1   y coords of the Ey/Hx sub-grid [m]
%
%                 Ex and Hy are co-located, as are Ey and Hx; MODE_POWER
%                 relies on exactly that pairing.  F.dx and F.dy are set by
%                 the callers after this function returns.  The coordinate
%                 vectors are used only by commented-out diagnostic plots,
%                 but keeping them is recommended.
%
% ---------------------------------------------------------------------------
% NOTES FOR ANYONE SUBSTITUTING A SOLVER
% ---------------------------------------------------------------------------
% HANDEDNESS.  The original solver returned fields in a LEFT-handed Cartesian
%   frame (it shared a code path with a cylindrical rho-phi-z solver), so its
%   [Ex Ey Ez Hx Hy Hz] correspond to [Ex -Ey Ez Hx -Hy Hz] in a standard
%   right-handed frame.  Everything downstream uses only |Ex|^2, |Ey|^2 and
%   |Ez|^2 at the sidewalls plus the guided power, so the sign convention
%   cancels and either handedness reproduces the paper.
%
% NORMALIZATION.  No particular field normalization is required.  The sidewall
%   field integrals are divided by the power returned by MODE_POWER, so any
%   consistent scaling of E and H cancels.
%
% HALF-PIXEL OFFSETS.  The callers index F.Ex, F.Ey and F.Ez with the SAME
%   (ix, iy) ranges when integrating along a sidewall, even though those three
%   arrays live on sub-grids offset by half a pixel and differ in size by one
%   element.  That is Tymon's original approximation, not an artifact of this
%   stub; it is second-order in the grid spacing and the published results
%   were produced with it.  A replacement solver should match the sub-grid
%   layout above rather than pre-interpolate everything onto a common grid,
%   or the stage-1 output will differ slightly from data/boundary_fields/.
%
% VALIDATION.  After plugging in a solver, run
%   examples/example4_validate_stage1.m
%   which regenerates part of the Si sweep and compares it cell-by-cell
%   against the shipped reference.
%
% See docs/modesolver_interface.md for the longer discussion.

error('waveguide_scattering:noModeSolver', ...
      ['No mode solver is configured.\n' ...
       'Stage 1 requires a 2D full-vector mode solver, which is not shipped ' ...
       'with this project.\nImplement modesolve_cross_section.m (and ' ...
       'mode_power.m) to run stage 1, or use the precomputed stage-1 data ' ...
       'in data/boundary_fields/ and start from stage 2.\n' ...
       'See docs/modesolver_interface.md.']);

end
