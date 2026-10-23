function []=scat_hic_integral_parallel_merge(serial_part_file, shard_files, out_file)
% Merge parallel stage-2 shards (plus an optional serial checkpoint) into the
% standard '<filename>_integrated' workspace file, indistinguishable to
% downstream scat_hic_interpolate_* scripts from a serial run's output.
%
%   serial_part_file : path to a serial checkpoint containing completed rows
%                      (pass '' if the whole sweep was computed by shards)
%   shard_files      : cell array of '<filename>_integrated_shard<id>' paths
%   out_file         : output path (the standard '<filename>_integrated' name)
%
% Coverage is verified: every filled cell (data(ia,ib).a > 0) must have a
% nonzero PaynePower_TE (strictly positive when computed).
% GNW / Claude, 2026-08-19.

CELL_VARS = {'PaynePower_TE','PowerHIC_TE_Y_cos_1','PowerHIC_TE_Y_cos_2', ...
             'PowerHIC_TE_Z_cos_1','PowerHIC_TE_Z_cos_2','PowerHIC_TE_YZ_cos_1', ...
             'PowerHIC_TE_YZ_cos_2','PowerHIC_TE_norm','PowerHIC_TE_comb', ...
             'PowerHIC_TE_test','PowerHIC_TE_factor','PowerHIC_TE_total', ...
             'PaynePower_TM','PowerHIC_TM_X_cos_1','PowerHIC_TM_X_cos_2', ...
             'PowerHIC_TM_Z_sin_1','PowerHIC_TM_Z_sin_2','PowerHIC_TM_norm', ...
             'PowerHIC_TM_factor','PowerHIC_TM_total', ...
             'matrix_a','matrix_b','matrix_phi_squared_TM','matrix_phi_squared_TE'};

if ~isempty(serial_part_file)
  S = load(serial_part_file);
else
  S = load(shard_files{1});   % seed with first shard (has data, dxy, constants)
  for v = CELL_VARS
    if isfield(S, v{1}), S.(v{1})(:) = 0; end
  end
end

nrow_total = size(S.data, 1);
ncol_total = size(S.data, 2);

for s = 1:numel(shard_files)
  H = load(shard_files{s});
  rows = H.ia_list(:).';
  for v = CELL_VARS
    name = v{1};
    if ~isfield(H, name)
      error('shard %s missing variable %s', shard_files{s}, name);
    end
    if ~isfield(S, name), S.(name) = zeros(nrow_total, ncol_total); end
    % pad target to full size if needed
    [r1, c1] = size(S.(name));
    S.(name)(max(r1,nrow_total), max(c1,ncol_total)) = S.(name)(r1, c1)*0 + S.(name)(r1, c1); %#ok<NASGU>
    [hr, hc] = size(H.(name));
    for ia = rows
      if ia <= hr
        S.(name)(ia, 1:hc) = H.(name)(ia, :);
      end
    end
  end
end

% coverage check
missing = 0;
for ia = 1:nrow_total
  for ib = 1:ncol_total
    d = S.data(ia,ib);
    if isfield(d,'a') && ~isempty(d.a) && any(d.a(:)>0)
      if S.PaynePower_TE(ia,ib) == 0
        fprintf('MISSING cell (%d,%d)\n', ia, ib);
        missing = missing + 1;
      end
    end
  end
end
if missing > 0
  error('merge incomplete: %d filled cells lack integrated values', missing);
end

S.savefilename = out_file;
save(out_file, '-struct', 'S');
fprintf('merged %d shards -> %s (coverage complete)\n', numel(shard_files), out_file);

end
