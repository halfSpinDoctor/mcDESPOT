% SCRIPT Test Parameter Sweep of Tau


% Simulation Parameters
SNR = Inf;
N   = 1; % Number of noise realizations

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

fract = 0:.001:.30;
% tau   = 0.075:0.0005:0.250;
tau   = .100;

% Calibration Maps here
omega = 20.00;  % Off-resonance,  B0
fam   = 1.0;    % Flip Angle Map, B1

for jj = 1:length(fract)

  % Generate Noise realizations
  for ii = 1:N
    % -- Generate Signal Curves Using CPU --
    [data_spgr(ii,:) data_ssfp_0(ii,:)  ] = mcdespot_sim_v2([1000 1 1/t1_f 1/t2_f 1./t1_m 1./t2_m fract(jj) tau omega*2*pi], alpha_spgr, alpha_ssfp, tr_spgr, tr_ssfp, SNR);
    [x               data_ssfp_180(ii,:)] = mcdespot_sim_v2([1000 1 1/t1_f 1/t2_f 1./t1_m 1./t2_m fract(jj) tau (omega*tr_ssfp*2*pi + pi)/tr_ssfp], alpha_spgr, alpha_ssfp, tr_spgr, tr_ssfp, SNR);
    
    [c_data_spgr c_data_ssfp_0 c_data_ssfp_180] = cpMCDESPOT_signal(t1_m, t1_f, t2_m, t2_f, fract(jj), tau, omega, alpha_spgr', alpha_ssfp', tr_spgr, tr_ssfp);
    
    plot(alpha_ssfp, data_ssfp_180, '-', alpha_ssfp, c_data_ssfp_180*(mean(data_ssfp_180(:))./mean(c_data_ssfp_180(:))), 'o');
    title(['MWF = ' num2str(fract(jj))]);
    pause(.01);
    
    M(ii) = getframe;
    
  end

  % [fv rnrm] = cpMCDESPOT_model_fit(data_spgr, data_ssfp_0, data_ssfp_180, alpha_spgr, alpha_ssfp, tr_spgr, tr_ssfp, fam*ones([N 1]), zeros([N 1]), mode, debug);
  

  
  % fv_tau_mean(jj,:) = mean(fv);
  % fv_tau_std(jj,:)  = std(fv);
  % fv = [];
end