% EXAMPLE 2 -- Run the model itself on a few waveguides.     (~30 s per cell)
%
% Example 1 read a precomputed answer out of a map.  This one actually runs
% the stage-2 physics: it takes the boundary-field summary of a waveguide
% cross-section and evaluates the radiated-power integrals of the paper
% (adaptive Lobatto quadrature over the far-field radiation profile), then
% checks the result against the shipped stage-2 output.
%
% This is the example to read if you want to see what the model consumes and
% what it produces.  No mode solver needed.

here = fileparts(mfilename('fullpath'));
root = fileparts(here);
addpath(fullfile(root, 'src', 'waveguide_scattering'));
addpath(fullfile(root, 'src', 'waveguide_scattering', 'parallel'));
warning('off', 'Octave:possible-matlab-short-circuit-operator');
warning('off', 'Octave:function-name-clash');

ROW    = 28;   % row 28 of the Si sweep is the a = 195 nm row, nearest the
               % 200-nm-high waveguide the paper uses as its worked example.
NCELLS = 4;    % how many widths along that row to compute (the full row is
               % 19 cells; raise this if you want the whole row).

% ---- what goes in ---------------------------------------------------------
S = load(fullfile(root, 'data', 'boundary_fields', 'boundary_fields_Si_v3.mat'));
d = S.data(ROW, 1);

fprintf('Stage-2 input for one waveguide cross-section:\n');
fprintf('  geometry    height a = %.1f nm, width b = %.1f nm\n', d.P.t*1e9, d.P.wd*1e9);
fprintf('  materials   n_core = %.2f, n_clad = %.2f, lambda = %.0f nm\n', ...
        d.P.nco, d.P.ncl, d.P.l0*1e9);
fprintf('  TE-like     beta  = %.6g rad/m  (n_eff = %.4f)\n', ...
        real(d.TE.beta), real(d.TE.beta)*d.P.l0/(2*pi));
fprintf('              Phi^2 = %.6g   (sidewall field intensity / guided power)\n', d.TE.Phi_squared);
fprintf('              R_YZ  = %.6g   (|Ez|^2 / |Ex|^2 on the sidewall)\n', d.TE.R_YZ);
fprintf('  TM-like     beta  = %.6g rad/m  (n_eff = %.4f)\n', ...
        real(d.TM.beta), real(d.TM.beta)*d.P.l0/(2*pi));
fprintf('              Phi^2 = %.6g\n', d.TM.Phi_squared);
fprintf('              R_XZ  = %.6g\n', d.TM.R_XZ);
fprintf('\nThat is the whole interface: a dozen numbers per cross-section, all\n');
fprintf('of them reductions of the modal field onto the two rough sidewalls.\n');
fprintf('See docs/data_format.md.\n\n');

% ---- stage 2 on a truncated copy of that row ------------------------------
% scat_hic_integral_final_v3 integrates an entire sweep.  The parallel worker
% is the same per-cell code with the row loop exposed, so it can do one row;
% and the row loop stops at the first cell whose .a is not positive, so
% blanking cell NCELLS+1 in a COPY of the data shortens the row for this demo.
rundir = fullfile(root, 'output', 'example2');
if ~exist(rundir, 'dir'), mkdir(rundir); end

data = S.data;
if NCELLS < size(data, 2)
  data(ROW, NCELLS + 1).a = 0;
end
save('-v7', fullfile(rundir, 'boundary_fields_Si_v3.mat'), 'data');

old = pwd();
cd(rundir);
fprintf('Integrating %d cells of row %d (~30 s each) ...\n\n', NCELLS, ROW);
t0 = tic;
scat_hic_integral_parallel_worker('boundary_fields_Si_v3', ROW, 'example2');
elapsed = toc(t0);
cd(old);

C = load(fullfile(rundir, 'boundary_fields_Si_v3_integrated_shardexample2'));
R = load(fullfile(root, 'data', 'integrated', 'boundary_fields_Si_v3_integrated.mat'));

% ---- what comes out -------------------------------------------------------
fprintf('\nDone in %.0f s.\n\n', elapsed);
fprintf('%-8s %-8s | %-14s %-14s | %-11s\n', ...
        'a (nm)', 'b (nm)', 'TE (dB/cm)/nm2', 'TM (dB/cm)/nm2', 'max rel err');
fprintf('%s\n', repmat('-', 1, 70));
worst = 0;
for ib = 1:NCELLS
  te = C.PowerHIC_TE_total(ROW, ib);
  tm = C.PowerHIC_TM_total(ROW, ib);
  ref = [R.PowerHIC_TE_total(ROW, ib), R.PowerHIC_TM_total(ROW, ib)];
  err = max(abs([te, tm] - ref) ./ abs(ref));
  worst = max(worst, err);
  fprintf('%-8.1f %-8.1f | %-14.5f %-14.5f | %.2e\n', ...
          C.matrix_a(ROW, ib)*1e9, C.matrix_b(ROW, ib)*1e9, te, tm, err);
end

fprintf('\nColumns 3-4 are the Fig. 12 quantity: multiply by the roughness\n');
fprintf('variance in nm^2 to get the propagation loss in dB/cm.\n');
fprintf('Column 5 compares against data/integrated/.  The quadrature is\n');
fprintf('deterministic, so agreement should be at floating-point level.\n');
if worst < 1e-9
  fprintf('Worst disagreement %.2e -- reproduced.\n', worst);
else
  fprintf('*** Worst disagreement %.2e is larger than expected. ***\n', worst);
end
