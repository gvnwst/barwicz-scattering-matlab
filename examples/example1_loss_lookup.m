% EXAMPLE 1 -- What the model predicts, in physical units.    (~1 minute)
%
% Runs stage 3 on the shipped Si/SiO2 results and reads a scattering loss out
% of the resulting map for one specific waveguide, then checks it against the
% number quoted in the caption of Fig. 12 of the paper:
%
%   "if one expects sigma^2 = 10 nm^2, the expected loss in the TE-like mode
%    of a 200-nm-high and 400-nm-wide Si waveguide will be around 7 dB/cm."
%
% No mode solver needed -- this starts from data/integrated/.

here = fileparts(mfilename('fullpath'));
root = fileparts(here);
addpath(fullfile(root, 'src', 'waveguide_scattering'));
warning('off', 'Octave:possible-matlab-short-circuit-operator');
warning('off', 'Octave:function-name-clash');
set(0, 'DefaultFigureVisible', 'off');

% ---- the waveguide and roughness we are asking about ----------------------
height_nm   = 200;    % core height  (a in the code, y in the paper)
width_nm    = 400;    % core width   (b in the code, x in the paper)
sigma_sq_nm2 = 10;    % roughness variance [nm^2]
%   Correlation length is fixed at Lc = 50 nm: it is baked into the shipped
%   stage-2 results.  Changing it means rerunning stage 2 (see example 3).

% ---- stage 3: coarse results -> interpolated loss map ---------------------
% scat_hic_interpolate_Si is a script-like function that loads a fixed
% filename from the current directory, so give it a scratch directory
% containing that file under the name it expects.
rundir = fullfile(root, 'output', 'example1');
if ~exist(rundir, 'dir'), mkdir(rundir); end
copyfile(fullfile(root, 'data', 'integrated', 'boundary_fields_Si_v3_integrated.mat'), ...
         fullfile(rundir, 'boundary_fields_Si_v3_integrated.mat'));

old = pwd();
cd(rundir);
fprintf('Running stage 3 (interpolation) for Si/SiO2 ...\n');
scat_hic_interpolate_Si;
cd(old);

M = load(fullfile(rundir, 'scat_hic_interpolate_Si'));

% ---- read the map at the requested cross-section --------------------------
nm = 1e-9;
loss_per_nm2 = interp2(M.mesh_b_i, M.mesh_a_i, M.mesh_PowerHIC_TE_total_i, ...
                       width_nm * nm, height_nm * nm, 'linear');
eff_per_nm2  = interp2(M.mesh_b_i, M.mesh_a_i, M.mesh_PowerHIC_TE_factor_i, ...
                       width_nm * nm, height_nm * nm, 'linear');

loss_dB_per_cm = loss_per_nm2 * sigma_sq_nm2;

fprintf('\n');
fprintf('Si/SiO2 waveguide, lambda = 1550 nm, exponential roughness, Lc = 50 nm\n');
fprintf('  core             : %g nm high x %g nm wide\n', height_nm, width_nm);
fprintf('  polarization     : TE-like (fundamental)\n');
fprintf('  roughness        : sigma^2 = %g nm^2\n', sigma_sq_nm2);
fprintf('\n');
fprintf('  Fig. 12 quantity : %.4f (dB/cm)/nm^2\n', loss_per_nm2);
fprintf('  Fig. 11 quantity : %.4f  (3-D loss / 2-D Payne-Lacey loss)\n', eff_per_nm2);
fprintf('  -> SCATTERING LOSS = %.2f dB/cm\n', loss_dB_per_cm);
fprintf('\n');
fprintf('  paper, Fig. 12 caption : "around 7 dB/cm"\n');
if abs(loss_dB_per_cm - 7) < 1
  fprintf('  AGREES with the published value.\n');
else
  fprintf('  *** does NOT match the published value -- investigate. ***\n');
end
fprintf('\n');
fprintf('The Fig. 11 quantity below 1 is the central result of the paper: for\n');
fprintf('this cross-section a 2-D analysis overestimates the loss %.1fx.\n', 1/eff_per_nm2);
