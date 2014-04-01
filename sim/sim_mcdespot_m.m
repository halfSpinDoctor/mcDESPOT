% FUNCTION [s_spgr s_ssfp_0 s_ssfp_180] = sim_mcdespot_m(fv, alpha, tr_spgr, tr_ssfp, snr)
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
%       v1.0 - initial version, fixed exp() --> expm() for matrix exponential (Sept-2009)
%       v4.0 Fixed Documentation, fixed 2*pi*TR hack for off-resonance

function [s_spgr s_ssfp_0 s_ssfp_180] = sim_mcdespot_m(fv, alpha, tr_spgr, tr_ssfp, snr)

% Convert flip angles to radians
alpha = alpha .* pi/180;

% Convert X to parameter values
r1_m    = 1./fv(1);
r1_f    = 1./fv(2);

r2_m    = 1./fv(3);
r2_f    = 1./fv(4);

MWF     = fv(5);
tau_f   = fv(6);

omega   = fv(7)*tr_ssfp*2*pi;

pd_spgr = fv(8);
pd_ssfp = fv(9);

k_sf = MWF./(tau_f*(1-MWF));
k_fs = 1./tau_f;
  
M_spgr = [];
  
for jj = 1:length(alpha)
  % -- Mz of SPGR (Eq [1] in Deoni 2008 MRM) --
    
  % *** Note, first term (1,1) of A_spgr has been corrected (quited as -1/(t1_f-k_fs) in paper)
  A_spgr = [-r1_m-k_fs    k_sf  
             k_fs        -r1_f-k_sf];

  disp(['Rank of A_SPGR:   ' num2str(rank(A_spgr))]);
  disp(['Det  of A_SPGR:   ' num2str(det( A_spgr))]);
  disp(['Cond of A_SPGR:   ' num2str(cond(A_spgr))]);
  disp(['Det  of exp(A*TR) ' num2str(det(expm(A_spgr*tr_spgr)))]);
  disp(['Cond of exp(A*TR) ' num2str(cond(expm(A_spgr*tr_spgr)))]);
  
  M0 = pd_spgr * [MWF (1-MWF)]';
  
  % Mss - Steady state longitudinal magnitization, times a factor sin(alpha) gives the 
  % observed Mxy signal (absorbs T2* and reciever gain into PD Term)
  M_spgr_full = (eye(2) - expm(A_spgr*tr_spgr) * cos(alpha(jj))) \ (eye(2) - expm(A_spgr.*tr_spgr)) * M0 * sin(alpha(jj)); 
  M_spgr(jj)  = M_spgr_full(1) + M_spgr_full(2); %#ok<AGROW>
  
end
  
  
% -- Mz of SSFP Eq [4] Deoni 2008) --
  
for jj = 1:length(alpha)

  % == SSFP-0 ==

  w = (omega + 0)/tr_ssfp;

  A_ssfp = [-r2_m-k_fs   k_sf       -w           0          0            0
            k_fs        -r2_f-k_sf   0          -w          0            0
            w            0          -r2_m-k_fs   k_sf       0            0
            0            w           k_fs       -r2_f-k_sf  0            0
            0            0           0           0         -r1_m-k_fs    k_sf
            0            0           0           0          k_fs        -r1_f-k_sf];

  A_ssfp;
  disp(['Rank of A_SSFP:   ' num2str(rank(A_ssfp))]);
  disp(['Det  of A_SSFP:   ' num2str(det( A_ssfp))]);
  disp(['Cond of A_SSFP:   ' num2str(cond(A_ssfp))]);
  disp(['Det  of exp(A*TR) ' num2str(det(expm(A_ssfp*tr_ssfp)))]);
  disp(['Cond of exp(A*TR) ' num2str(cond(expm(A_ssfp*tr_ssfp)))]);

  % Rotation Matrix
  R = [1    0    0                 0               0               0
       0    1    0                 0               0               0
       0    0    cos(alpha(jj))    0              -sin(alpha(jj))  0
       0    0    0                 cos(alpha(jj))  0              -sin(alpha(jj))
       0    0    sin(alpha(jj))    0               cos(alpha(jj))  0
       0    0    0                 sin(alpha(jj))  0               cos(alpha(jj))];
   
  C = pd_ssfp .* [0 0 0 0 MWF (1-MWF)]';

  % Full magnitization vector
  % M = [Mx_f Mx_s My_f My_s Mz_f Mz_s]  (Use Eq[13] from MRM 2007 Investigating the Effect of Exchange and Multi...)
  M_ssfp_full = (eye(6) - expm(A_ssfp*tr_ssfp)*R) \ (eye(6) - expm(A_ssfp*tr_ssfp))*C;
 
  % Measured magnitization signal
  M_ssfp_0(jj) = abs(M_ssfp_full(1) + M_ssfp_full(2) + 1i*M_ssfp_full(3) + 1i*M_ssfp_full(4));    %#ok<AGROW>


  % == SSFP-180 ==
  w = (omega + pi)/tr_ssfp;
  
    A_ssfp = [-r2_m-k_fs   k_sf       -w           0          0            0
              k_fs        -r2_f-k_sf   0          -w          0            0
              w            0          -r2_m-k_fs   k_sf       0            0
              0            w           k_fs       -r2_f-k_sf  0            0
              0            0           0           0         -r1_m-k_fs    k_sf
              0            0           0           0          k_fs        -r1_f-k_sf];
            
            

  % Rotation Matrix
  R = [1    0    0                 0               0               0
       0    1    0                 0               0               0
       0    0    cos(alpha(jj))    0              -sin(alpha(jj))  0
       0    0    0                 cos(alpha(jj))  0              -sin(alpha(jj))
       0    0    sin(alpha(jj))    0               cos(alpha(jj))  0
       0    0    0                 sin(alpha(jj))  0               cos(alpha(jj))];
   
  C = pd_ssfp .* [0 0 0 0 MWF (1-MWF)]';

  % Full magnitization vector
  % M = [Mx_f Mx_s My_f My_s Mz_f Mz_s]  (Use Eq[13] from MRM 2007 Investigating the Effect of Exchange and Multi...)
  M_ssfp_full = (eye(6) - expm(A_ssfp*tr_ssfp)*R) \ (eye(6) - expm(A_ssfp*tr_ssfp))*C;
 
  % Measured magnitization signal (180)
  M_ssfp_180(jj) = abs(M_ssfp_full(1) + M_ssfp_full(2) + 1i*M_ssfp_full(3) + 1i*M_ssfp_full(4)); %#ok<AGROW>
  % M_ssfp_180(jj) = norm([M_ssfp_full(1) M_ssfp_full(3)]) + norm([M_ssfp_full(2) M_ssfp_full(4)]); %#ok<AGROW>

end
  
% Add Noise Gaussian Random Noise
max_signal = max([pd_spgr pd_ssfp]);
noise = max_signal/snr;
  
% DEBUG: Print out noise level
% disp(['Noise: ' num2str(noise)]);
  
s_spgr     = M_spgr     + noise*randn;
s_ssfp_0   = M_ssfp_0   + noise*randn;
s_ssfp_180 = M_ssfp_180 + noise*randn;
