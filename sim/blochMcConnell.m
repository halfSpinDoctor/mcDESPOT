% FUNCTION [dM] = blochMcConnell(M, fv, b1)
%
% Inputs:
%    M  - Two-pool magnetization Vector [Mx_m Mx_f My_m My_f Mz_m Mz_f]
%    fv - two-pool system parameters [t1_m t1_f t2_m t2_f MWF tau omega pd_spgr pd_ssfp] 
%         system paramters T1,2, Tau in s.  Omega off-resonance in Hz
%    b1 - radiofrequency excitation pulse envlope, in Tesla
%
% Outputs:
%   dM - Two-pool differential magnetization vector [dMx_m dMx_f dMy_m dMy_f dMz_m dMz_f]
%
% Samuel A. Hurley
% University of Wisconsin
% v1.0 25-Oct-2010
%
% Changelog:
%   1.0 - Initial version, based off of Lenz et al. Mag Res Mat Phys 2010
%

function dM = blochMcConnell(M, fv, b1)

% Constants
GAMMA = 267.513e6; % rad/s/T   1H Gyromagnetic Ratio

% Compute omega_0 and omega_1 in rad/s
w_0 = fv(7)*2*pi;  % 2pi to convert off-resonance from Hz to Rad/s
w_1 = GAMMA * b1;

% Pull Out NMR Parameters
% Convert X to parameter values
T1_m    = fv(1);
T1_f    = fv(2);

T2_m    = fv(3);
T2_f    = fv(4);

MWF     = fv(5);
tau     = fv(6);

M0      = fv(8);

% Compute K_mf and K_fm from MWF and tau
Kmf     = 1.00/tau;
Kfm     = Kmf * MWF/(1-MWF);

% Mag vector
Mx_m = M(1);
Mx_f = M(2);
My_m = M(3);
My_f = M(4);
Mz_m = M(5);
Mz_f = M(6);

% == Differential ==
dM(1) =   w_0*My_m   - Mx_m/T2_m  - Kmf*Mx_m + Kfm*Mx_f;
dM(2) =   w_0*My_f   - Mx_f/T2_f  - Kfm*Mx_f + Kmf*Mx_m;

dM(3) =  -w_0*Mx_m   - My_m/T2_m  - Kmf*My_m + Kfm*My_f  + w_1*Mz_m;
dM(4) =  -w_0*Mx_f   - My_f/T2_f  - Kfm*My_f + Kmf*My_m  + w_1*Mz_f;

dM(5) =  (M0*MWF     - Mz_m)/T1_m - Kmf*Mz_m + Kfm*Mz_f  - w_1*My_m;
dM(6) =  (M0*(1-MWF) - Mz_f)/T1_f - Kfm*Mz_f + Kmf*Mz_m  - w_1*My_f;

dM = dM';

return;