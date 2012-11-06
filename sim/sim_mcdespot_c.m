% FUNCTION [s_spgr s_ssfp_0 s_ssfp_180] = sim_mcdespot_c2(fv, alpha, tr_spgr, tr_ssfp, snr)
%
% Simulation for SAH Version of C-Code (Re-written with John Ollinger's modifications)
%
% T1m = T1 fast/myelin   T1f = T1 slow/"free"
%
% Inputs:
%    fv               - Parameter Vector [t1_m t1_f t2_m t2_f MWF tau omega pd_spgr pd_ssfp]
%                       t1, t2, tau in seconds, 0 < MWF < 1, omega in Hz, pd in a.u.
%    alpha            - flip angles [degrees]
%    tr_spgr, tr_ssfp - repitition times
%    snr              - signal to noise ratio rel to maximum proton density term (spgr or ssfp)
%
% Outputs:
%    s_spgr
%    s_ssfp_0
%    s_ssfp_180
%
% Samuel A. Hurley
% University of Wisconsin
% v4.0 11-Nov-2010
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
  % For Sam's c-code, only 3 outputs & FV is now a matrix, not a set of vectors
  % [resSPGR resSSFP_0 resSSFP_180] = cpMCDESPOT_residuals_SAH(fv(1:7), 0, 0, 0, alpha(ii), alpha(ii), tr_spgr, tr_ssfp, 1); 
  
  % Extract signals from the residual
  s_spgr(ii)     = sqrt(cpMCDESPOT_residuals_SAH(fv(1:6)', fv(7), -1,   0, alpha(ii), tr_spgr, 1));
  s_ssfp_0(ii)   = sqrt(cpMCDESPOT_residuals_SAH(fv(1:6)', fv(7),  0,   0, alpha(ii), tr_ssfp, 1));
  s_ssfp_0(ii)   = sqrt(cpMCDESPOT_residuals_SAH(fv(1:6)', fv(7),  180, 0, alpha(ii), tr_ssfp, 1));
end

% Add Noise Gaussian Random Noise
max_signal = 1;
noise = max_signal/snr;
  
% DEBUG: Print out noise level
% disp(['Noise: ' num2str(noise)]);
  
s_spgr     = s_spgr     + noise*randn;
s_ssfp_0   = s_ssfp_0   + noise*randn;
s_ssfp_180 = s_ssfp_180 + noise*randn;
