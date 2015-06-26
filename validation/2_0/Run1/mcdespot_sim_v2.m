% FUNCTION [data_spgr data_ssfp] = mcdespot_sim(x, flip_spgr, flip_ssfp, tr_spgr, tr_ssfp, noise)
%
% Function to generate simulated mcDESPOT data for testing minimization
% routine with gaussian random noise
%
% Inputs:
%        x - [pd_spgr, pd_ssfp, r1_s, r2_s, r1_f, r2_f, f_f, tau_f]
%            (Simulated Parameter Values, r1&2 in seconds)
%
%        tr_spr, tr_ssfp - TR of sequence (in s)
%        fmap - off-resonance effect (in Hz)
%        noise - SNR of zero-mean gaussian random noise (rel. to max signal observed)
%
% Outputs:
%
%        data_spgr, data_ssfp - simulated scanner data (arbitrary units)
%
% Pouria Mossahebi
% Samuel A. Hurley
% University of Wisconsin
% v2.1 2-Feb-2010
%
% Revision History:
%     v2.0 - implemented diagonalization of A in order to
%            improve performance of matrix exponential.
% 
%          - changed input from individual parameters to a single x vector
%            for simplified function calling  [18-Sept-2009]
%
%     v2.1  -changed to rician noise, fixed random noise (before same offset
%            was added to all points in the curve!)
%
%           -changed variable fmap to omega
%

function [data_spgr data_ssfp] = mcdespot_sim_v2(x, alpha_spgr, alpha_ssfp, tr_spgr, tr_ssfp, snr)

% Convert flip angles to radians
alpha_spgr = alpha_spgr .* pi/180;
alpha_ssfp = alpha_ssfp .* pi/180;

% Convert X to parameter values
pd_spgr = x(1);
pd_ssfp = x(2);
r1_s    = x(3);
r2_s    = x(4);
r1_f    = x(5);
r2_f    = x(6);
f_f     = x(7);
tau_f   = x(8);
omega   = x(9);

k_sf = f_f./(tau_f*(1-f_f));
k_fs = 1./tau_f;

M_spgr = [];

for jj = 1:length(alpha_spgr)
  % -- Mz of SPGR (Eq [1] in Deoni 2008 MRM) --

  % *** Note, first term (1,1) of A_spgr has been modified from -1/(t1_f-k_fs)
  A_spgr = [-r1_f-k_fs    k_sf  
             k_fs        -r1_s-k_sf];


  M0 = pd_spgr * [f_f (1-f_f)]';

  % Mss - Steady state longitudinal magnitization, times a factor sin(alpha) gives the observed Mxy signal (ignoring T2* decay)
  
  % Diagonalize A to improve matrix exponential performance
  [U L] = eig(A_spgr);
  
  % Take e^A*TR for diagonal terms only
  Exp_L_tr_spgr = diag(exp(diag(L*tr_spgr)));
  
  % Full magnitization vector
  M_spgr_full = U * inv(eye(2) - Exp_L_tr_spgr * cos(alpha_spgr(jj)))* (eye(2) - Exp_L_tr_spgr) * inv(U) * M0 * sin(alpha_spgr(jj));
  M_spgr(jj)  = M_spgr_full(1) + M_spgr_full(2);

end


% -- Mz of SSFP Eq [4] Deoni 2008) --
M_ssfp = [];

for jj = 1:length(alpha_ssfp)
  A_ssfp = [-r2_f-k_fs   k_sf        omega       0          0            0
            k_fs        -r2_s-k_sf   0           omega      0            0
           -omega        0          -r2_f-k_fs   k_sf       0            0
            0           -omega       k_fs       -r2_s-k_sf  0            0
            0            0           0           0         -r1_f-k_fs    k_sf
            0            0           0           0          k_fs        -r1_s-k_sf];

  % Rotation Matrix
  R = [1    0    0                    0                   0                   0
       0    1    0                    0                   0                   0
       0    0    cos(alpha_ssfp(jj))  0                   sin(alpha_ssfp(jj)) 0
       0    0    0                    cos(alpha_ssfp(jj)) 0                   sin(alpha_ssfp(jj))
       0    0    -sin(alpha_ssfp(jj)) 0                   cos(alpha_ssfp(jj)) 0
       0    0    0                   -sin(alpha_ssfp(jj)) 0                   cos(alpha_ssfp(jj))];

  C = pd_ssfp .* [0 0 0 0 f_f.*r1_f (1-f_f)*r1_s];

  % Diagonalize A to improve matrix exponential performance
  [U L] = eig(A_ssfp);    
  Inv_U = inv(U);
  Exp_L_tr_ssfp = diag(exp(diag(L*tr_ssfp)));
  
  % Full Magnitization vector
  % M = [Mx_f Mx_s My_f My_s Mz_f Mz_s]  (Use Eq[13] from MRM 2007 Investigating the Effect of Exchange and Multi...)
  M_ssfp_full = U * inv(eye(6)-Exp_L_tr_ssfp*Inv_U*R*U)*(Exp_L_tr_ssfp - eye(6)) * inv(L)* Inv_U * C';

  % Measured magnitization signal
  M_ssfp(jj) = abs(M_ssfp_full(1) + M_ssfp_full(2) + i*M_ssfp_full(3) + i*M_ssfp_full(4));

end

% Add Noise Gaussian Random Noise
max_signal = max([M_spgr(:); M_ssfp(:)]);
noise = max_signal/snr;

% DEBUG: Print out noise level
% disp(['Noise: ' num2str(noise)]);

data_spgr = abs(M_spgr + noise*randn(size(M_spgr)) + 1i*noise*randn(size(M_spgr)));
data_ssfp = abs(M_ssfp + noise*randn(size(M_ssfp)) + 1i*noise*randn(size(M_ssfp)));

end
