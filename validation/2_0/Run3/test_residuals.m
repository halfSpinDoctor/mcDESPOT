% SCRIPT test_mcdespot
% Script to test the performance of mcDESPOT fitting on changing one parameter
% at a time, using Monte Carlo Noise testing

% Samuel A. Hurley
% Pouria Mossahebi
% University of Wisconsin
% v1.0 2-Feb

tic;

% Movie output options
MOV = 0;  % 1 for AVI movie

% Simulation Parameters
SNR = Inf;
N   = 50; % Number of noise realizations

% Fitting Mode
mode  = 0;  % Full fit 180+0+B0 fitting
debug = 0;

% Scan Parameters
tr_spgr = .0065;
tr_ssfp = .0050;

% Use flip angles from our actual scans
alpha_spgr = [3 4 5 6 7 9 13 18];
alpha_ssfp = [12 16 21 27 33 40 51 68];


% ----   #6  ----

% T1_m (T1 myelin) = T1_f (in paper - T1 FAST)
% T1_f (T1 FREE)   = T1_s (in paper - T1 SLOW)

t1_m = 0.465;
t1_f = 1.07;

t2_m = .026;
t2_f = .117;

fract      = 0.20;
mean_tau   = 0.075:.001:0.250;
tau        = 0.075:.0005:0.250;

% Calibration Maps here
omega = 20.00;  % Off-resonance,  B0
fam   = 1.0;    % Flip Angle Map, B1

if MOV == 1
  mov=avifile('SSFP_180_Tau_Res','compression','none'); 
end


for ii = 1:length(mean_tau)


  for jj = 1:length(tau)
    
    x(jj,:) = [t1_m t1_f t2_m t2_f fract tau(jj) omega];
    
  end

[data_spgr data_ssfp_0  ] = mcdespot_sim_v2([1000 1 1/t1_f 1/t2_f 1./t1_m 1./t2_m fract mean_tau(ii) omega(1)*2*pi], alpha_spgr, alpha_ssfp, tr_spgr, tr_ssfp, SNR);
[~,        data_ssfp_180] = mcdespot_sim_v2([1000 1 1/t1_f 1/t2_f 1./t1_m 1./t2_m fract mean_tau(ii) (omega(1)*tr_ssfp*2*pi + pi)/tr_ssfp], alpha_spgr, alpha_ssfp, tr_spgr, tr_ssfp, SNR);


[resSPGR resSSFP_0 resSSFP_180 medSSFP_0 medSSFP_180] = cpMCDESPOT_residuals(x(:,1), x(:,2), x(:,3), x(:,4), x(:,5), x(:,6), x(:,7)*tr_ssfp*2*pi, data_spgr'/mean(data_spgr(:)), data_ssfp_0'/mean(data_ssfp_0(:)), data_ssfp_180'/mean(data_ssfp_180(:)), alpha_spgr', alpha_ssfp', tr_spgr, tr_ssfp);

plot(tau, resSSFP_0+resSSFP_180+resSPGR);
title(['Overall Residual Actual Tau: ' num2str(mean_tau(ii))]);
xlim([0.075 0.250]);
%ylim([0 1200]);

xlabel 'Tau [s]'
ylabel 'Overall Residual [a.u.]'

drawnow;

if MOV == 1
  f2=getframe(gcf); % gets the gcf  
  mov=addframe(mov,f2); % adds the frame into mov  
end


end

if MOV == 1
  mov=close(mov); % closes the mov  
end


toc;