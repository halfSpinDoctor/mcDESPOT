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

%%

clc
tic;

% Simulation Parameters
SNR = Inf;

% Scan Parameters
tr_spgr = .0050;             % s
trs     = .001:.001:.020;    % increments of SSFP TR

alpha       = 1:.1:68; % degrees
alpha_bloch = 2:6:50;

% -- Default Tissue Parameters (Defined by Deoni in paper).
% T1_m (T1 myelin) = T1_f (in paper - T1 FAST)
% T1_f (T1 FREE)   = T1_s (in paper - T1 SLOW)
pd_spgr = 1;
pd_ssfp = 1;

t1_m  = 0.465;        % s
t1_f  = 1.07;         % s

t2_m  = 0.026;        % s
t2_f  = 0.117;        % s

MWF   = 0.30;         % 20%
tau   = 0.200;        % s
omega = 020;            % Range of off-resonance freq


%% Loop to test how curve moves /w different TR
for ii = 1:length(trs)
  fv = [t1_m t1_f t2_m t2_f MWF tau omega pd_spgr pd_ssfp];
  
  tr_ssfp = trs(ii);
  
  disp('--------------');
  tic
  [m2_spgr m2_ssfp_0 m2_ssfp_180]    = sim_mcdespot_m( fv, alpha, tr_spgr, tr_ssfp, SNR);       % M Code
  toc
  [c_spgr c_ssfp_0 c_ssfp_180]       = sim_mcdespot_c( fv, alpha, tr_spgr, tr_ssfp, SNR);       % Seans C Code
  toc
  [c2_spgr c2_ssfp_0 c2_ssfp_180]    = sim_mcdespot_c2(fv, alpha, tr_spgr, tr_ssfp, SNR);       % My C Code
  toc
  [s_spgr s_ssfp_0 s_ssfp_180]       = sim_mcdespot_bloch_2(fv, alpha_bloch, tr_spgr, tr_ssfp, 1e-8); % My Bloch Simulator
  toc;
  
  
  
  %% Plot SD C-Code vs Bloch
  plot(alpha, c_ssfp_180, alpha, c_ssfp_0, alpha_bloch, s_ssfp_180, 'rx', alpha_bloch, s_ssfp_0, 'bx')
  legend('SCLD-SSFP-180', 'SCLD-SSFP-0', 'Bloch SSFP-180', 'Bloch SSFP-0');
  xlabel 'Flip Angle [degrees]'
  ylabel 'MR Signal [a.u.]'
  ylim([0 .16]);
  title(['mcDESPOT ' num2str(MWF*100) '% MWF, TR_S_S_F_P = ' num2str(tr_ssfp)]);
  
  savefig(['SD_TR' num2str(tr_ssfp*1000, '%02.0f') '.tif']);
  
  close all;
  
  %% Plot SAH C-Code vs Bloch
  plot(alpha, c2_ssfp_180, alpha, c2_ssfp_0, alpha_bloch, s_ssfp_180, 'rx', alpha_bloch, s_ssfp_0, 'bx')
  legend('SAH-SSFP-180', 'SAH-SSFP-0', 'Bloch SSFP-180', 'Bloch SSFP-0');
  xlabel 'Flip Angle [degrees]'
  ylabel 'MR Signal [a.u.]'
  ylim([0 .16]);
  title(['mcDESPOT ' num2str(MWF*100) '% MWF, TR_S_S_F_P = ' num2str(tr_ssfp)]);

  savefig(['SAH_TR' num2str(tr_ssfp*1000, '%02.0f') '.tif']);
  
  close all;
  
  
end


