% SCRIPT test_mcdespot_trmsim
% Script to test Deoni's code vs. Sam's new re-write of C code
% /w John Ollinger's matrix exponential implementation.
%
% Samuel A. Hurley
% University of Wisconsin
% v4.0 11-Nov-2010
% v4.1 15-Nov-2010
% v4.2 01-Feb-2011
%      07-Feb-2011
% v5.0 12-May-2011  Example to show S.D. at ISMRM
% v6.0 07-Mar-2014  For TobiasWood Paper

%%

clc
tic;

% Simulation Parameters
SNR = Inf;

% Scan Parameters
alpha   = 20; % degrees

tr_spgr = .005;
tr_ssfp = 0.0060;
% tau     = 0.0032; % t_{rf}
tau     = 0;
t1      = 1.100;
t2      = 0.060;
pd      = 1;

omega   = -165:1:165;

s_ssfp = [];


%% Loop over off-resonance Omega to plot frequency response function
for ii = 1:length(omega)
  
  fv = [t1 t2 omega(ii) pd]; % FV for single-pool (Bloch) equations
  
  disp('--------------');
  tic
  [x s_ssfp(ii) x m_ssfp] = sim_mcdespot_bloch_single(fv, alpha, tr_spgr, tr_ssfp, tau, 1e-4); % My Bloch Simulator
  toc
  
end

save _run_03
  
%   %% Plot SD C-Code vs Bloch
%   plot(alpha, c_ssfp_180, alpha, c_ssfp_0, alpha_bloch, s_ssfp_180, 'rx', alpha_bloch, s_ssfp_0, 'bx')
%   legend('SCLD-SSFP-180', 'SCLD-SSFP-0', 'Bloch SSFP-180', 'Bloch SSFP-0');
%   xlabel 'Flip Angle [degrees]'
%   ylabel 'MR Signal [a.u.]'
%   ylim([0 .16]);
%   title(['mcDESPOT ' num2str(MWF*100) '% MWF, TR_S_S_F_P = ' num2str(tr_ssfp)]);
%   
%   savefig(['SD_TR' num2str(tr_ssfp*1000, '%02.0f') '.tif']);
%   
%   close all;
%   
%   %% Plot SAH C-Code vs Bloch
%   plot(alpha, c2_ssfp_180, alpha, c2_ssfp_0, alpha_bloch, s_ssfp_180, 'rx', alpha_bloch, s_ssfp_0, 'bx')
%   legend('SAH-SSFP-180', 'SAH-SSFP-0', 'Bloch SSFP-180', 'Bloch SSFP-0');
%   xlabel 'Flip Angle [degrees]'
%   ylabel 'MR Signal [a.u.]'
%   ylim([0 .16]);
%   title(['mcDESPOT ' num2str(MWF*100) '% MWF, TR_S_S_F_P = ' num2str(tr_ssfp)]);
% 
%   savefig(['SAH_TR' num2str(tr_ssfp*1000, '%02.0f') '.tif']);
%   
%   close all;
%   
%   



