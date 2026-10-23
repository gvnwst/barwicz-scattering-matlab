function [P, SzTM, SzTE] = mode_power(F, m)
% MODE_POWER  Guided power of a mode, split into its two polarization terms
%             (NOT IMPLEMENTED -- stage-1 plug-in point).
%
% Companion stub to MODESOLVE_CROSS_SECTION.  The original implementation was
% part of the same mode-solver package, which is not included here; see
% that file's header and docs/modesolver_interface.md.
%
% ---------------------------------------------------------------------------
% WHAT IT COMPUTES
% ---------------------------------------------------------------------------
% The z-directed guided power of mode m, as the integral of (E x H) . zhat
% over the cross-section, with the two cross-product terms kept separate:
%
%       SzTE = -Ex(:,:,m) .* Hy(:,:,m)          (MX-1) x (MY-2)
%       SzTM = +Ey(:,:,m) .* Hx(:,:,m)          (MX-2) x (MY-1)
%       P    = sum( [SzTE(:); SzTM(:)] ) * F.dx * F.dy
%
% Three things about this differ from the textbook expression and are
% deliberate -- a replacement must reproduce them to stay consistent with
% data/boundary_fields/:
%
%   1. NO COMPLEX CONJUGATE on H, and no factor of 1/2.  The original code
%      used the unconjugated product.  For the lossless bound modes of the
%      paper the fields can be taken real up to a common phase, so this
%      agrees with (1/2)Re{E x H*} up to an overall real factor -- which
%      cancels, since P only ever appears as a denominator.
%   2. THE SIGN of the SzTE term is negative, compensating the left-handed
%      frame the original solver worked in (see MODESOLVE_CROSS_SECTION).
%   3. THE TWO TERMS LIVE ON DIFFERENT SUB-GRIDS and are deliberately not
%      interpolated onto a common one.  Ex and Hy are co-located, and so are
%      Ey and Hx, so each product is a valid elementwise multiply; the two
%      resulting arrays simply have different shapes.  Summing both and
%      multiplying by one dx*dy is the original's midpoint-rule integral.
%
% ---------------------------------------------------------------------------
% INPUTS / OUTPUTS
% ---------------------------------------------------------------------------
%   F      field struct from MODESOLVE_CROSS_SECTION, carrying F.dx, F.dy.
%   m      mode index to evaluate (the callers pass 1, then 2).
%
%   P      scalar guided power of mode m, in whatever units follow from the
%          field normalization.  Used only as the denominator of
%          TE.Phi_squared / TM.Phi_squared, so its absolute scale is
%          arbitrary but must be consistent with the E fields in F.
%   SzTM   the +Ey.*Hx array described above.
%   SzTE   the -Ex.*Hy array described above.
%          The caller forms sum(sum(abs(SzTE))) and sum(sum(abs(SzTM))) for
%          both modes and compares them to decide which mode is TE-like and
%          which is TM-like; only their relative magnitude matters there.
%
% CONSISTENCY REQUIREMENT
%   P must be the power of the SAME field arrays the caller then integrates
%   along the sidewalls.  If your solver rescales or interpolates fields,
%   compute P from the final arrays, or Phi_squared picks up a spurious
%   scale factor and the predicted loss is wrong by that factor.
%
% See docs/modesolver_interface.md.

error('waveguide_scattering:noModeSolver', ...
      ['No mode-power routine is configured.\n' ...
       'Implement mode_power.m together with modesolve_cross_section.m to ' ...
       'run stage 1,\nor use the precomputed stage-1 data in ' ...
       'data/boundary_fields/ and start from stage 2.\n' ...
       'See docs/modesolver_interface.md.']);

end
