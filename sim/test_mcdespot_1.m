% SCRIPT test_mcdespot
% Script to test the accuracy of the mcDESPOT objective function vs.
% Full Bloch Differential Equation Simulation
%
% Samuel A. Hurley
% University of Wisconsin
% v4.0 25-Oct-2010

tic;

% Simulation Parameters
SNR = Inf;

% Scan Parameters
tr_spgr = .0065;
tr_ssfp = .0050;
alpha   = .5:.5:90;

% -- Default Tissue Parameters (Defined by Deoni in paper).
% T1_m (T1 myelin) = T1_f (in paper - T1 FAST)
% T1_f (T1 FREE)   = T1_s (in paper - T1 SLOW)
pd_spgr = 1;
pd_ssfp = 1;

t1_m = 0.465;
t1_f = 1.07;

t2_m = .026;
t2_f = .117;

MWF  = .20;
tau  = 0.200; % Lit value quoted in Lenz & Scheffler
omega = 0:1:200;    % Perfectly on-resonance

for ii = 1:length(omega)
  fv = [t1_m t1_f t2_m t2_f MWF tau omega(ii) pd_spgr pd_ssfp];
  
  [m_spgr m_ssfp_0 m_ssfp_180] = sim_mcdespot_m(fv, alpha, tr_spgr, tr_ssfp, SNR);
  [c_spgr c_ssfp_0 c_ssfp_180] = sim_mcdespot_c(fv, alpha, tr_spgr, tr_ssfp, SNR);
  
  subplot(1,2,1);
  plot(alpha, m_spgr, alpha, m_ssfp_0, alpha, m_ssfp_180, alpha, c_spgr, alpha, c_ssfp_0, alpha, c_ssfp_180);
  legend('SPGR M', 'SSFP-0 M', 'SSFP-180 M', 'SPGR C', 'SSFP-0 C', 'SSFP-180 C');
  xlabel 'Flip Angle [degrees]'
  ylabel 'MR Signal [a.u.]'
  title 'Signal Curves: C and M Code';
  
  subplot(1,2,2);
  plot(alpha, m_spgr-c_spgr, alpha, m_ssfp_0-c_ssfp_0, alpha, m_ssfp_180-c_ssfp_180, alpha, repmat(eps, [1 length(alpha)]), ':k', alpha, repmat(-eps, [1 length(alpha)]), ':k');
  legend('SPGR', 'SSFP-0', 'SSFP-180', '+/- Machine Eps');
  xlabel 'Flip Angle [degrees]'
  ylabel 'MR Signal [a.u.]'
  title 'Signal Curves: M-C Code';
  
  pause(.250);
end

% saveas(gcf, 'mcDESPOT_C_vs_M_Test.tif');

