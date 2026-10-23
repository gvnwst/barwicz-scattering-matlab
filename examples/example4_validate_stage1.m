% EXAMPLE 4 -- Validate a mode solver you have plugged in.   (~2 min per row)
%
% Stage 1 is the only part of the pipeline that needs a mode solver, and no
% mode solver ships with this project (see docs/modesolver_interface.md).
% Once you have implemented modesolve_cross_section.m and mode_power.m, run
% this to check your implementation: it regenerates a few rows of the Si
% sweep and compares every stored quantity against data/boundary_fields/,
% which was produced with the solver Tymon Barwicz originally used.
%
% If you have not plugged in a solver yet, this script says so and stops.

here = fileparts(mfilename('fullpath'));
root = fileparts(here);
addpath(fullfile(root, 'src', 'waveguide_scattering'));
warning('off', 'Octave:possible-matlab-short-circuit-operator');
warning('off', 'Octave:function-name-clash');
set(0, 'DefaultFigureVisible', 'off');

ROWS = 1:3;   % rows of the Si sweep to regenerate and check

% ---- is a solver actually installed? --------------------------------------
try
  modesolve_cross_section(struct(), 1, 1, struct(), 1, [0 0 0 0], [0 0]);
  have_solver = true;
catch err
  have_solver = ~strcmp(err.identifier, 'waveguide_scattering:noModeSolver');
  if ~have_solver
    fprintf('No mode solver is plugged in, so stage 1 cannot run.\n\n');
    fprintf('%s\n\n', err.message);
    fprintf('Everything else in this repository works without one:\n');
    fprintf('  example1_loss_lookup.m            read a loss out of the maps\n');
    fprintf('  example2_single_cell.m            run the scattering model\n');
    fprintf('  example3_reproduce_fig11_fig12.m  reproduce the paper figures\n');
    return
  end
end

% ---- regenerate part of the Si sweep --------------------------------------
% boundary_fields_Si_v3.m sweeps all 32 rows and saves after every cell, so
% run it and stop once the rows of interest are on disk, or edit its max_a /
% min_a bounds for a shorter sweep.
rundir = fullfile(root, 'output', 'example4');
if ~exist(rundir, 'dir'), mkdir(rundir); end

old = pwd();
cd(rundir);
fprintf('Running stage 1 with your mode solver ...\n');
boundary_fields_Si_v3;
cd(old);

% ---- compare against the reference ----------------------------------------
G = load(fullfile(rundir, 'boundary_fields_Si_v3.mat'));
R = load(fullfile(root, 'data', 'boundary_fields', 'boundary_fields_Si_v3.mat'));

fields = {'TE.beta', 'TE.R_YZ', 'TE.Phi_squared', ...
          'TM.beta', 'TM.R_XZ', 'TM.Phi_squared'};

fprintf('\n%-8s %-8s %-16s %-12s\n', 'a (nm)', 'b (nm)', 'quantity', 'rel. diff');
fprintf('%s\n', repmat('-', 1, 52));
worst = 0;
for ia = ROWS
  for ib = 1:size(R.data, 2)
    if isempty(R.data(ia, ib).P), break; end
    for f = 1:numel(fields)
      parts = strsplit(fields{f}, '.');
      got = G.data(ia, ib).(parts{1}).(parts{2});
      ref = R.data(ia, ib).(parts{1}).(parts{2});
      rel = abs(got - ref) / abs(ref);
      worst = max(worst, rel);
      if rel > 1e-6
        fprintf('%-8.1f %-8.1f %-16s %.3e\n', R.data(ia, ib).a*1e9, ...
                R.data(ia, ib).b*1e9, fields{f}, rel);
      end
    end
  end
end

fprintf('\nWorst relative difference over rows %s: %.3e\n', mat2str(ROWS), worst);
if worst < 1e-6
  fprintf('Your solver reproduces the reference stage-1 data.\n');
  fprintf('(Only differences above 1e-6 are listed above, so the table is empty.)\n');
else
  fprintf(['Differences above 1e-6 are listed above. Some are expected: near the\n' ...
           'square cross-sections (a == b) the two modes are degenerate and the\n' ...
           'TE/TM assignment is arbitrary, so those cells can differ by O(1).\n' ...
           'See docs/validation.md.\n']);
end
