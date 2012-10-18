% SCRIPT test_mcdespot_trmsim
% Script to test Deoni's code vs. Sam's new re-write of C code
% /w John Ollinger's matrix exponential implementation.
%
% Samuel A. Hurley
% University of Wisconsin
% v4.0 11-Nov-2010
% v4.1 15-Nov-2010


clc
tic;

% Simulation Parameters
SNR = Inf;

% Scan Parameters
tr_spgr = .0050;      % s
%tr_ssfp = .0100;      % s
trs = .005;
alpha   = 2:6:46; % degrees

% -- Default Tissue Parameters (Defined by Deoni in paper).
% T1_m (T1 myelin) = T1_f (in paper - T1 FAST)
% T1_f (T1 FREE)   = T1_s (in paper - T1 SLOW)
pd_spgr = 1;
pd_ssfp = 1;

t1_m  = 0.465;        % s
t1_f  = 1.07;         % s

t2_m  = .026;         % s
t2_f  = .117;         % s

MWF   = .40;
tau   = 0.800;        % Lit value quoted in Lenz & Scheffler
omega = 0;            % Range of off-resonance freq


% Loop to test how curve moves /w off-resonance
for ii = 1:length(trs)
  fv = [t1_m t1_f t2_m t2_f MWF tau omega pd_spgr pd_ssfp];
  tr_ssfp = trs(ii);
  
  disp('--------------');
  tic
  [m_spgr m_ssfp_0 m_ssfp_180]    = sim_mcdespot_m( fv, alpha, tr_spgr, tr_ssfp, SNR);
  
  fv(5) = .5;
  
  [m2_spgr m2_ssfp_0 m2_ssfp_180]    = sim_mcdespot_m( fv, alpha, tr_spgr, tr_ssfp, SNR);
  toc
  [c_spgr c_ssfp_0 c_ssfp_180]    = sim_mcdespot_c( fv, alpha, tr_spgr, tr_ssfp, SNR);
  toc
  [c2_spgr c2_ssfp_0 c2_ssfp_180] = sim_mcdespot_c2(fv, alpha, tr_spgr, tr_ssfp, SNR);
  toc
  
  [s_spgr s_ssfp_0 s_ssfp_180] = sim_mcdespot_bloch_2(fv, alpha, tr_spgr, tr_ssfp, 1e-8);
  toc;
  


%% Plot M-code vs Bloch
subplot(1,2,1);
plot(alpha, m_ssfp_180, alpha, m_ssfp_0, alpha, m_spgr, alpha, m2_ssfp_180, alpha, m2_ssfp_0, alpha, m2_spgr)
% legend('M-SSFP-180', 'M-SSFP-0', 'SPGR');
% legend('M-SSFP-180', 'Bloch-SSFP-180', 'M-SSFP-0', 'Bloch-SSFP-0', 'M-SPGR', 'Bloch-SPGR');
xlabel 'Flip Angle [degrees]'
ylabel 'MR Signal [a.u.]'
ylim([0 .16]);
title(['mcDESPOT, TE=TR, 40% vs 50% MWF']);
subplot(1,2,2);
plot(alpha, m_ssfp_180 - m2_ssfp_180, alpha, m_ssfp_0 - m2_ssfp_0, alpha, m_spgr - m2_spgr)
% legend('M-SSFP-180', 'M-SSFP-0', 'SPGR');
% legend('M-SSFP-180', 'Bloch-SSFP-180', 'M-SSFP-0', 'Bloch-SSFP-0', 'M-SPGR', 'Bloch-SPGR');
xlabel 'Flip Angle [degrees]'
ylabel 'MR Signal [a.u.]'
ylim([-.1 .1]);
title(['Differences (40% vs 50% MWF)']);
legend('SSFP-180', 'SSFP-0', 'SPGR');

saveas(gcf, 'MWF_Bias.tif')


end


