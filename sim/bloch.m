% FUNCTION [dM] = bloch(M, fv, b1)
%
% Inputs:
%    M  - Single-pool magnetization Vector [Mx My Mz]
%    fv - single-pool system parameters [t1 t2 omega pd] 
%         system paramters T1,2 in s.  Omega off-resonance in Hz (1/s)
%    b1 - radiofrequency excitation pulse envlope, in Tesla
%
% Outputs:
%   dM - Single-pool differential magnetization vector [dMx dMy dMz]
%
% Samuel A. Hurley
% University of Wisconsin
% v1.0 31-Jan-2011
%
% Changelog:
%   1.0 - Initial version, based off of blochMcConnell.m
%

function dM = bloch(M, fv, b1)

% Constants
GAMMA = 267.513e6; % rad/s/T   1H Gyromagnetic Ratio

% Compute omega_0 and omega_1 in rad/s
w_0 = fv(3)*2*pi;  % 2pi to convert off-resonance from Hz to Rad/s
w_1 = GAMMA * b1;

% Pull Out NMR Parameters
% Convert X to parameter values
T1    = fv(1);
T2    = fv(2);

M0    = fv(4);

% Mag vector
Mx = M(1);
My = M(2);
Mz = M(3);

% == Differential ==
dM(1) =  w_0*My - Mx/T2;

dM(2) = -w_0*Mx - My/T2 + w_1*Mz_m;

dM(3) = (M0 - Mz)/T1    - w_1*My_m;

dM = dM';

return;