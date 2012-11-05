% FUNCTION [s_spgr s_ssfp_0 s_ssfp_180] = sim_mcdespot(fv, alpha, tr_spgr, tr_ssfp, snr)
%
% Simulation for SD Version of C-Code (Re-written with John Ollinger's modifications)
%
% T1m = T1 fast/myelin   T1f = T1 slow/"free"
%
% Inputs:
%    fv               - Parameter Vector [t1_m t1_f t2_m t2_f MWF tau omega pd_spgr pd_ssfp]
%                       t1, t2, tau in seconds, 0 < MWF < 1, omega in Hz, pd in a.u.
%    alpha            - flip angles
%    tr_spgr, tr_ssfp - repitition times
%    snr              - signal to noise ratio rel to maximum proton density term (spgr or ssfp)
%
% Outputs:
%    s_spgr
%    s_ssfp_0
%    s_ssfp_180
%
% Samuel A. Hurley
% Pouria Mossahebi
% University of Wisconsin
% v4.0 22-Oct-2010
%
% Changelog:
%       v1.0 Initial Version, Jun-2010
%       v4.0 Fixed Documentation, added gaussian noise

function [s_spgr s_ssfp_0 s_ssfp_180] = sim_mcdespot_c(fv, alpha, tr_spgr, tr_ssfp, snr)

% Preallocate some stuff
s_spgr     = zeros([1 length(alpha)]);
s_ssfp_0   = zeros([1 length(alpha)]);
s_ssfp_180 = zeros([1 length(alpha)]);

% Loop over flip angles
for ii = 1:length(alpha)
  [resSPGR resSSFP_0 resSSFP_180 x x] = cpMCDESPOT_residuals_SD(fv(1), fv(2), fv(3), fv(4), fv(5), fv(6), fv(7)*tr_ssfp*2*pi, fv(8), fv(9), 0, 0, 0, alpha(ii), alpha(ii), tr_spgr, tr_ssfp); %#ok<NASGU>
  % Extract signals from the residual
  s_spgr(ii) = sqrt(resSPGR);
  s_ssfp_0(ii) = sqrt(resSSFP_0);
  s_ssfp_180(ii) = sqrt(resSSFP_180);
end

% Add Noise Gaussian Random Noise
max_signal = max(fv(8:9));
noise = max_signal/snr;
  
% DEBUG: Print out noise level
% disp(['Noise: ' num2str(noise)]);
  
s_spgr     = s_spgr     + noise*randn;
s_ssfp_0   = s_ssfp_0   + noise*randn;
s_ssfp_180 = s_ssfp_180 + noise*randn;
