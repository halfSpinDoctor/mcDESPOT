% FUNCTION [s_spgr s_ssfp_0 s_ssfp_180] = sim_mcdespot(fv, omega, alpha, tr_spgr, tr_ssfp, snr)
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
% v1.1 13-Sep-2015
%
% Changelog:
%       v1.0 Initial Version, based off of sim_mcdespot_c2
%            Alows for multiple FV vectors to be passed to a single sim_mcdespot call
%       v1.1 Updated to pass numThreads argument for updated C-code

function [s_spgr s_ssfp_0 s_ssfp_180] = sim_mcdespot(fv, omega, alpha, tr_spgr, tr_ssfp, snr)

% Preallocate some stuff
s_spgr     = zeros([size(fv,1) length(alpha)]);
s_ssfp_0   = zeros([size(fv,1) length(alpha)]);
s_ssfp_180 = zeros([size(fv,1) length(alpha)]);

% Loop over flip angles
for ii = 1:length(alpha)
  % For Sam's c-code, only 3 outputs & FV is now a matrix, not a set of vectors
  resSPGR = cpMCDESPOT_residuals_SAH(fv', omega,   -1, 0, alpha(ii), tr_spgr,1 , 1);
  resSSFP_0 = cpMCDESPOT_residuals_SAH(fv', omega,    0, 0, alpha(ii), tr_ssfp,1 , 1);
  resSSFP_180 = cpMCDESPOT_residuals_SAH(fv', omega,  180, 0, alpha(ii), tr_ssfp,1 , 1);
  
  % Extract signals from the residual
  s_spgr(:, ii) = sqrt(resSPGR);
  s_ssfp_0(:, ii) = sqrt(resSSFP_0);
  s_ssfp_180(:, ii) = sqrt(resSSFP_180);
end

% Add Noise Gaussian Random Noise
max_signal = 1;
noise = max_signal/snr;
  
% DEBUG: Print out noise level
% disp(['Noise: ' num2str(noise)]);
  
s_spgr     = s_spgr     + noise*randn;
s_ssfp_0   = s_ssfp_0   + noise*randn;
s_ssfp_180 = s_ssfp_180 + noise*randn;
