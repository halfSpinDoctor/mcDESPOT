% SCRIPT test_mcdespot_c_vs_c2
% Script to test Deoni's code vs. Sam's new re-write of C code
% /w John Ollinger's matrix exponential implementation.
%
% Samuel A. Hurley
% University of Wisconsin
% v4.0 11-Nov-2010

clc
tic;

% Simulation Parameters
SNR = Inf;

% Scan Parameters
tr_spgr = .0050;      % s
tr_ssfp = .0050;      % s
alpha   = 1:.5:89.5; % degrees

% -- Default Tissue Parameters (Defined by Deoni in paper).
% T1_m (T1 myelin) = T1_f (in paper - T1 FAST)
% T1_f (T1 FREE)   = T1_s (in paper - T1 SLOW)
pd_spgr = 1;
pd_ssfp = 1;

t1_m  = 0.465;        % s
t1_f  = 1.07;         % s

t2_m  = .026;         % s
t2_f  = .117;         % s

MWF   = .20;
tau   = 0.200;        % Lit value quoted in Lenz & Scheffler
omega = 0;            % Range of off-resonance freq


% Loop to test how curve moves /w off-resonance
for ii = 1:length(t2_f)
  fv = [t1_m t1_f t2_m t2_f MWF tau omega pd_spgr pd_ssfp];
  
  disp('--------------');
  tic
  [m_spgr m_ssfp_0 m_ssfp_180]    = sim_mcdespot_m( fv, alpha, tr_spgr, tr_ssfp, SNR);
  toc
  [c_spgr c_ssfp_0 c_ssfp_180]    = sim_mcdespot_c( fv, alpha, tr_spgr, tr_ssfp, SNR);
  toc
  [c2_spgr c2_ssfp_0 c2_ssfp_180] = sim_mcdespot_c2(fv, alpha, tr_spgr, tr_ssfp, SNR);
  toc
  
%   [b_spgr b_ssfp_0 b_ssfp_180]    = sim_mcdespot_bloch(fv, alpha, tr_spgr, tr_ssfp, 2.4e-14);
%   toc
  
  
  %% Plot comparison of results
  subplot(1,2,1);
  plot(alpha, c_spgr-m_spgr, alpha, c_ssfp_0-m_ssfp_0, alpha, c_ssfp_180-m_ssfp_180, alpha, repmat(eps, [1 length(alpha)]), ':k', alpha, repmat(-eps, [1 length(alpha)]), ':k');
  legend('SPGR', 'SSFP-0', 'SSFP-180', '+/- Machine Eps');
  xlabel 'Flip Angle [degrees]'
  ylabel 'MR Signal [a.u.]'
  % ylim([-2.5e-3 0]);
  title 'Signal Curves: C Code vs M Code';
  
  subplot(1,2,2);
  plot(alpha, c2_spgr-m_spgr, alpha, c2_ssfp_0-m_ssfp_0, alpha, c2_ssfp_180-m_ssfp_180, alpha, repmat(eps, [1 length(alpha)]), ':k', alpha, repmat(-eps, [1 length(alpha)]), ':k');
  legend('SPGR', 'SSFP-0', 'SSFP-180', '+/- Machine Eps');
  xlabel 'Flip Angle [degrees]'
  ylabel 'MR Signal [a.u.]'
  % ylim([-2.5e-3 0]);
  title 'Signal Curves: New C Code vs M Code';


end

% Save comparison plot
%saveas(gcf, 'M_vs_C_vs_C2_Comparison.tif');

