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
tr_spgr = .005;
tr_ssfp = .001;
alpha   = [10];

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
% omega = 0:1:200;    % Perfectly on-resonance

omega = 80; % Hz
% NOTE: Changing the off-resonance seems to change the condition number of A, but not
% the condition number of expm(A). Need to investigate further.

for ii = 1:length(omega)
  fv = [t1_m t1_f t2_m t2_f MWF tau omega(ii) pd_spgr pd_ssfp];
  
  % [m_spgr m_ssfp_0 m_ssfp_180] = sim_mcdespot_m(fv, alpha, tr_spgr, tr_ssfp, SNR);
  [m_spgr m_ssfp_0 m_ssfp_180] = sim_mcdespot_m(fv, alpha, tr_spgr, tr_ssfp, 1e-6);
  
  pause(.250);

end

% savefig('/Users/samuel/Dropbox/Fig.png');
% 
% NOTIFY_EMAIL = '2628943063@email.uscc.net'; 
% send_mail_message(NOTIFY_EMAIL, 'ProcDone', ['Processing Complete at ' datetime()]);
