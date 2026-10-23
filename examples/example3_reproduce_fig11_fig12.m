% EXAMPLE 3 -- Reproduce Figs. 11 and 12 of the paper.         (~2-5 minutes)
%
% Runs stage 3 for all three material systems on the shipped stage-2 results
% and writes the interpolated loss maps to output/example3/.  Those maps are
% the six panels of Fig. 11 (radiation efficiency) and the six panels of
% Fig. 12 (loss normalized to roughness variance).
%
% Rendering: this script attempts the contour plots in Octave/MATLAB itself.
% The published-quality composites in figures/ were rendered separately by
% examples/plot_fig11_fig12.py (matplotlib), which reads the mesh files this
% script writes; use that if your Octave build cannot print headless.
%
% No mode solver needed.

here = fileparts(mfilename('fullpath'));
root = fileparts(here);
addpath(fullfile(root, 'src', 'waveguide_scattering'));
warning('off', 'Octave:possible-matlab-short-circuit-operator');
warning('off', 'Octave:function-name-clash');
set(0, 'DefaultFigureVisible', 'off');

rundir = fullfile(root, 'output', 'example3');
if ~exist(rundir, 'dir'), mkdir(rundir); end

% name, integrated dataset, stage-3 script, axis limits in nm, panel letters
systems = { ...
  {'SiON',  'boundary_fields_SiON_v3_e_integrated',  'scat_hic_interpolate_SiON',  [750 4000], 'ab', 'SiON/SiO2 (n_{core}=1.50)'}, ...
  {'Si3N4', 'boundary_fields_Si3N4_v3_d_integrated', 'scat_hic_interpolate_Si3N4', [200 1000], 'cd', 'Si_3N_4/SiO2 (n_{core}=2.00)'}, ...
  {'Si',    'boundary_fields_Si_v3_integrated',      'scat_hic_interpolate_Si',    [120  600], 'ef', 'Si/SiO2 (n_{core}=3.50)'} };

meshes = struct();

for k = 1:numel(systems)
  name = systems{k}{1};  dataset = systems{k}{2};  script = systems{k}{3};

  fprintf('[%d/3] %-6s  stage 3 ...\n', k, name);
  copyfile(fullfile(root, 'data', 'integrated', [dataset '.mat']), ...
           fullfile(rundir, [dataset '.mat']));

  old = pwd();
  cd(rundir);
  feval(script);
  cd(old);
  close all   % the stage-3 scripts leave their own diagnostic figures open

  M = load(fullfile(rundir, script));
  % flat variable names (e.g. Si_te_total) so the .mat is trivial to read
  % from MATLAB, Octave or scipy.io.loadmat alike
  meshes.([name '_a'])         = M.mesh_a_i;
  meshes.([name '_b'])         = M.mesh_b_i;
  meshes.([name '_te_total'])  = M.mesh_PowerHIC_TE_total_i;
  meshes.([name '_tm_total'])  = M.mesh_PowerHIC_TM_total_i;
  meshes.([name '_te_factor']) = M.mesh_PowerHIC_TE_factor_i;
  meshes.([name '_tm_factor']) = M.mesh_PowerHIC_TM_factor_i;
  clear M
end

meshfile = fullfile(rundir, 'stage3_meshes.mat');
save('-v7', meshfile, '-struct', 'meshes');
fprintf('\nInterpolated maps written to %s\n', meshfile);

% ---- render ---------------------------------------------------------------
% Octave's OpenGL (qt/fltk) toolkits frequently print black frames when no
% display is attached; the gnuplot toolkit does not. Prefer it if present.
if exist('graphics_toolkit', 'file')
  if any(strcmp('gnuplot', available_graphics_toolkits()))
    warning('off', 'Octave:gnuplot-graphics');
    graphics_toolkit('gnuplot');
  else
    fprintf(['\nNote: the gnuplot toolkit is unavailable. If the PNGs below ' ...
             'come out blank,\nrender them instead with ' ...
             '"python examples/plot_fig11_fig12.py".\n']);
  end
end

printed = true;
for fig = 1:2
  if fig == 1
    keys  = {'te_factor', 'tm_factor'};
    fname = 'fig11_reproduction.png';
  else
    keys  = {'te_total', 'tm_total'};
    fname = 'fig12_reproduction.png';
  end

  h = figure('Position', [0 0 1000 1200]);
  for k = 1:numel(systems)
    name = systems{k}{1};  lims = systems{k}{4};
    tags = systems{k}{5};  title_str = systems{k}{6};
    bb = meshes.([name '_b']) * 1e9;
    aa = meshes.([name '_a']) * 1e9;
    for c = 1:2
      subplot(3, 2, 2*(k-1) + c);
      z = meshes.([name '_' keys{c}]);
      % Line contours only, unlabelled: the maps are NaN outside the
      % single-mode region, which gnuplot's filled contours render badly, and
      % inline labels collide where the contours bunch. This is a quick look;
      % plot_fig11_fig12.py produces the filled, labelled version.
      contour(bb, aa, z, 12, 'k', 'LineWidth', 0.4);
      xlim(lims); ylim(lims);
      if c == 1, pol = 'TE-like'; else, pol = 'TM-like'; end
      title(sprintf('(%s) %s - %s', tags(c), title_str, pol), 'FontSize', 8);
      xlabel('Width of the core (nm)', 'FontSize', 8);
      ylabel('Height of the core (nm)', 'FontSize', 8);
    end
  end
  try
    print(h, fullfile(rundir, fname), '-dpng', '-r150');
    fprintf('wrote %s\n', fullfile(rundir, fname));
  catch err
    printed = false;
    fprintf('could not print %s: %s\n', fname, err.message);
  end
  close(h);
end

if ~printed
  fprintf('\nThis build could not print at all.\n');
end

fprintf('\nThose PNGs are a quick look: unlabelled line contours.\n');
fprintf('For the filled, labelled composites -- the ones in figures/ -- run\n');
fprintf('    python examples/plot_fig11_fig12.py\n');
fprintf('which reads %s.\n', meshfile);
fprintf('\nEither way, compare against Figs. 11 and 12 of the paper.\n');
