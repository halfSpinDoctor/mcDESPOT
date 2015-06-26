% SCRIPT test_mcdespot
% Script to test the performance of mcDESPOT fitting on changing one parameter
% at a time, using Monte Carlo Noise testing

% Samuel A. Hurley
% University of Wisconsin
% v3.0 10-Oct-2013

tic;

% Debug Plotting?
DEBUG = 0;

% Simulation Parameters
SNR   = 500; % Noise level
N     = 2000;  % Number of noise realizations

% Scan Parameters
tr_spgr = .0065;
tr_ssfp = .0050;

% Use flip angles from actual mcDESPOT v2 Protocol
alpha_spgr = [ 3  4  5  6  7  9 13 18];
alpha_ssfp = [12 16 21 27 33 40 51 68];

% ----   #1 - MWF Sweep, Hold Others Const  ----

% T1_m (T1 myelin) = T1_f (in paper - T1 FAST)
% T1_f (T1 FREE)   = T1_s (in paper - T1 SLOW)

t1_m = 0.465;
t1_f = 1.07;

t2_m = .026;
t2_f = .117;

mwf  = 0.00:0.02:0.30;
tau  = 0.180;

% % SD.MRM
% t1_m = 0.465;
% t1_f = 0.965;
% 
% t2_m = 0.012;
% t2_f = 0.090;
% 
% mwf  = 0.00:0.02:0.30;
% tau  = 0.125;

% Calibration Maps here
omega   = 40.00;  % Off-resonance,  B0
fam     = 1.0;    % Flip Angle Map, B1
pd_spgr = 1;
pd_ssfp = 1;

for jj = 1:length(mwf);
  % mcd Parameter Vector
  mcd_fv = [t1_m t1_f t2_m t2_f mwf(jj) tau omega pd_spgr pd_ssfp];

  % Generate Noise realizations
  
  % Preallocate
  data_spgr = zeros([N 8]);
  data_ssfp_0 = zeros([N 8]);
  data_ssfp_180 = zeros([N 8]);
  
  for ii = 1:N
    progressbar(ii/N);
    % -- Generate Signal Curves Using CPU --
     [data_spgr(ii,:) x x]                     = sim_mcdespot_c(mcd_fv, alpha_spgr, tr_spgr, tr_ssfp, SNR);
     [x data_ssfp_0(ii,:) data_ssfp_180(ii,:)] = sim_mcdespot_c(mcd_fv, alpha_ssfp, tr_spgr, tr_ssfp, SNR);
  end
  progressbar(1);

  [fv rnrm] = mcdespot_model_fit(data_spgr, data_ssfp_0, data_ssfp_180, alpha_spgr, alpha_ssfp, tr_spgr, tr_ssfp, fam*ones([N 1]), omega*ones([N 1]), repmat([t1_f t2_f pd_spgr pd_ssfp], [N 1]), DEBUG);
  
  fv_mwf_mean(jj,:) = mean(fv); %#ok<SAGROW>
  fv_mwf_std(jj,:)  = std(fv);  %#ok<SAGROW>
  
  fv = [];
  mcd_fv = [];
end


toc;

plot(mwf, fv_mwf_mean(:,5), '-k', mwf, fv_mwf_mean(:,5)-fv_mwf_std(:,5), 'k:', mwf, fv_mwf_mean(:,5)+fv_mwf_std(:,5), 'k:')
hold on

plot(mwf_x, fv_mwf_mean(2:15,5)'-sigma, 'r',  mwf_x, fv_mwf_mean(2:15,5)'+sigma, 'r');

% tmp = fv_mwf_mean(2:15,5);
% 
% errorbar(mwf_x(3:3:end), tmp(3:3:end), sigma_biased(3:3:end), 'b');

plot(mwf_x, fv_mwf_mean(2:15,5)'-sigma_biased_mean, 'b--',  mwf_x, fv_mwf_mean(2:15,5)'+sigma_biased_mean, 'b--');

for ii = 1:6
  sigma_param(ii)        = sqrt(var_min(ii,ii,10)) ./ fv_mwf_mean(11,ii);
  sigma_param_biased(ii) = sqrt(var_min_biased(ii,ii,10)) ./ fv_mwf_mean(11,ii);
  sigma_mc(ii)           = fv_mwf_std(11,ii) ./ fv_mwf_mean(11,ii);
  
end

data = [sigma_param; sigma_param_biased; sigma_mc];

bar(data', 'k')
set(gca,'YScale','log');
set(gca,'XTickLabel',['T1_M';'T1_F';'T2_M';'T2_F';'MWF ';'Tau '])


