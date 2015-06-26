% SCRIPT test_mcdespot
% Script to test the performance of mcDESPOT fitting on changing one parameter
% at a time, using Monte Carlo Noise testing

% Samuel A. Hurley
% Pouria Mossahebi
% University of Wisconsin
% v1.0 2-Feb

tic;

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

% ----   #1  ----


% T1_m (T1 myelin) = T1_f (in paper - T1 FAST)
% T1_f (T1 FREE)   = T1_s (in paper - T1 SLOW)

t1_m = 0.300:.001:.700;
t1_f = 1.07;

t2_m = .026;
t2_f = .117;

fract = 0.20;
tau   = 0.100;

% Calibration Maps here
omega = 20.00;  % Off-resonance,  B0
fam   = 1.0;    % Flip Angle Map, B1

for jj = 1:length(t1_m)

  % Generate Noise realizations
  for ii = 1:N
    % -- Generate Signal Curves Using CPU --
    [data_spgr(ii,:) data_ssfp_0(ii,:)  ] = mcdespot_sim_v2([1000 1 1/t1_f 1/t2_f 1./t1_m(jj) 1./t2_m fract tau omega*2*pi], alpha_spgr, alpha_ssfp, tr_spgr, tr_ssfp, SNR);
    [x               data_ssfp_180(ii,:)] = mcdespot_sim_v2([1000 1 1/t1_f 1/t2_f 1./t1_m(jj) 1./t2_m fract tau (omega*tr_ssfp*2*pi + pi)/tr_ssfp], alpha_spgr, alpha_ssfp, tr_spgr, tr_ssfp, SNR);
  end

  [fv rnrm] = cpMCDESPOT_model_fit(data_spgr, data_ssfp_0, data_ssfp_180, alpha_spgr, alpha_ssfp, tr_spgr, tr_ssfp, fam*ones([N 1]), zeros([N 1]), mode, debug);
  
  fv_t1m_mean(jj,:) = mean(fv);
  fv_t1m_std(jj,:)  = std(fv);
  fv = [];
end


% ----   #2  ----

% T1_m (T1 myelin) = T1_f (in paper - T1 FAST)
% T1_f (T1 FREE)   = T1_s (in paper - T1 SLOW)

t1_m = .465;
t1_f = .800:0.005:2.800;

t2_m = .026;
t2_f = .117;

fract = 0.20;
tau   = 0.100;

% Calibration Maps here
omega = 20.00;  % Off-resonance,  B0
fam   = 1.0;    % Flip Angle Map, B1

for jj = 1:length(t1_f)

  % Generate Noise realizations
  for ii = 1:N
    % -- Generate Signal Curves Using CPU --
    [data_spgr(ii,:) data_ssfp_0(ii,:)  ] = mcdespot_sim_v2([1000 1 1/t1_f(jj) 1/t2_f 1./t1_m 1./t2_m fract tau omega*2*pi], alpha_spgr, alpha_ssfp, tr_spgr, tr_ssfp, SNR);
    [x               data_ssfp_180(ii,:)] = mcdespot_sim_v2([1000 1 1/t1_f(jj) 1/t2_f 1./t1_m 1./t2_m fract tau (omega*tr_ssfp*2*pi + pi)/tr_ssfp], alpha_spgr, alpha_ssfp, tr_spgr, tr_ssfp, SNR);
  end

  [fv rnrm] = cpMCDESPOT_model_fit(data_spgr, data_ssfp_0, data_ssfp_180, alpha_spgr, alpha_ssfp, tr_spgr, tr_ssfp, fam*ones([N 1]), zeros([N 1]), mode, debug);
  

  
  fv_t1f_mean(jj,:) = mean(fv);
  fv_t1f_std(jj,:)  = std(fv);
  fv = [];
end


% ----   #3  ----


% T1_m (T1 myelin) = T1_f (in paper - T1 FAST)
% T1_f (T1 FREE)   = T1_s (in paper - T1 SLOW)

t1_m = .465;
t1_f = 1.07;

t2_m = .001:.0001:.040;
t2_f = .117;

fract = 0.20;
tau   = 0.100;

% Calibration Maps here
omega = 20.00;  % Off-resonance,  B0
fam   = 1.0;    % Flip Angle Map, B1

for jj = 1:length(t2_m)

  % Generate Noise realizations
  for ii = 1:N
    % -- Generate Signal Curves Using CPU --
    [data_spgr(ii,:) data_ssfp_0(ii,:)  ] = mcdespot_sim_v2([1000 1 1/t1_f 1/t2_f 1./t1_m 1./t2_m(jj) fract tau omega*2*pi], alpha_spgr, alpha_ssfp, tr_spgr, tr_ssfp, SNR);
    [x               data_ssfp_180(ii,:)] = mcdespot_sim_v2([1000 1 1/t1_f 1/t2_f 1./t1_m 1./t2_m(jj) fract tau (omega*tr_ssfp*2*pi + pi)/tr_ssfp], alpha_spgr, alpha_ssfp, tr_spgr, tr_ssfp, SNR);
  end

  [fv rnrm] = cpMCDESPOT_model_fit(data_spgr, data_ssfp_0, data_ssfp_180, alpha_spgr, alpha_ssfp, tr_spgr, tr_ssfp, fam*ones([N 1]), zeros([N 1]), mode, debug);
  

  
  fv_t2m_mean(jj,:) = mean(fv);
  fv_t2m_std(jj,:)  = std(fv);
  fv = [];
end

% ----   #4  ----


% T1_m (T1 myelin) = T1_f (in paper - T1 FAST)
% T1_f (T1 FREE)   = T1_s (in paper - T1 SLOW)

t1_m = 0.465;
t1_f = 1.07;

t2_m = .026;
t2_f = 0.05:0.0002:0.140;

fract = 0.20;
tau   = 0.100;

% Calibration Maps here
omega = 20.00;  % Off-resonance,  B0
fam   = 1.0;    % Flip Angle Map, B1

for jj = 1:length(t2_f)

  % Generate Noise realizations
  for ii = 1:N
    % -- Generate Signal Curves Using CPU --
    [data_spgr(ii,:) data_ssfp_0(ii,:)  ] = mcdespot_sim_v2([1000 1 1/t1_f 1/t2_f(jj) 1./t1_m 1./t2_m fract tau omega*2*pi], alpha_spgr, alpha_ssfp, tr_spgr, tr_ssfp, SNR);
    [x               data_ssfp_180(ii,:)] = mcdespot_sim_v2([1000 1 1/t1_f 1/t2_f(jj) 1./t1_m 1./t2_m fract tau (omega*tr_ssfp*2*pi + pi)/tr_ssfp], alpha_spgr, alpha_ssfp, tr_spgr, tr_ssfp, SNR);
  end

  [fv rnrm] = cpMCDESPOT_model_fit(data_spgr, data_ssfp_0, data_ssfp_180, alpha_spgr, alpha_ssfp, tr_spgr, tr_ssfp, fam*ones([N 1]), zeros([N 1]), mode, debug);

  
  fv_t2f_mean(jj,:) = mean(fv);
  fv_t2f_std(jj,:)  = std(fv);
  fv = [];
end

% ----   #5  ----


% T1_m (T1 myelin) = T1_f (in paper - T1 FAST)
% T1_f (T1 FREE)   = T1_s (in paper - T1 SLOW)

t1_m = 0.465;
t1_f = 1.07;

t2_m = .026;
t2_f = .117;

fract = 0.001:0.0005:0.250;
tau   = 0.100;

% Calibration Maps here
omega = 20.00;  % Off-resonance,  B0
fam   = 1.0;    % Flip Angle Map, B1

for jj = 1:length(fract)

  % Generate Noise realizations
  for ii = 1:N
    % -- Generate Signal Curves Using CPU --
    [data_spgr(ii,:) data_ssfp_0(ii,:)  ] = mcdespot_sim_v2([1000 1 1/t1_f 1/t2_f 1./t1_m 1./t2_m fract(jj) tau omega*2*pi], alpha_spgr, alpha_ssfp, tr_spgr, tr_ssfp, SNR);
    [x               data_ssfp_180(ii,:)] = mcdespot_sim_v2([1000 1 1/t1_f 1/t2_f 1./t1_m 1./t2_m fract(jj) tau (omega*tr_ssfp*2*pi + pi)/tr_ssfp], alpha_spgr, alpha_ssfp, tr_spgr, tr_ssfp, SNR);
  end

  [fv rnrm] = cpMCDESPOT_model_fit(data_spgr, data_ssfp_0, data_ssfp_180, alpha_spgr, alpha_ssfp, tr_spgr, tr_ssfp, fam*ones([N 1]), zeros([N 1]), mode, debug);
  

  
  fv_fract_mean(jj,:) = mean(fv);
  fv_fract_std(jj,:)  = std(fv);
  fv = [];
end

% ----   #6  ----


% T1_m (T1 myelin) = T1_f (in paper - T1 FAST)
% T1_f (T1 FREE)   = T1_s (in paper - T1 SLOW)

t1_m = 0.465;
t1_f = 1.07;

t2_m = .026;
t2_f = .117;

fract = 0.20;
tau   = 0.075:0.0005:0.250;

% Calibration Maps here
omega = 20.00;  % Off-resonance,  B0
fam   = 1.0;    % Flip Angle Map, B1

for jj = 1:length(tau)

  % Generate Noise realizations
  for ii = 1:N
    % -- Generate Signal Curves Using CPU --
    [data_spgr(ii,:) data_ssfp_0(ii,:)  ] = mcdespot_sim_v2([1000 1 1/t1_f 1/t2_f 1./t1_m 1./t2_m fract tau(jj) omega*2*pi], alpha_spgr, alpha_ssfp, tr_spgr, tr_ssfp, SNR);
    [x               data_ssfp_180(ii,:)] = mcdespot_sim_v2([1000 1 1/t1_f 1/t2_f 1./t1_m 1./t2_m fract tau(jj) (omega*tr_ssfp*2*pi + pi)/tr_ssfp], alpha_spgr, alpha_ssfp, tr_spgr, tr_ssfp, SNR);
  end

  [fv rnrm] = cpMCDESPOT_model_fit(data_spgr, data_ssfp_0, data_ssfp_180, alpha_spgr, alpha_ssfp, tr_spgr, tr_ssfp, fam*ones([N 1]), zeros([N 1]), mode, debug);
  

  
  fv_tau_mean(jj,:) = mean(fv);
  fv_tau_std(jj,:)  = std(fv);
  fv = [];
end

% ----   #7  ----


% T1_m (T1 myelin) = T1_f (in paper - T1 FAST)
% T1_f (T1 FREE)   = T1_s (in paper - T1 SLOW)

t1_m = 0.465;
t1_f = 1.07;

t2_m = .026;
t2_f = .117;

fract = 0.20;
tau   = 0.100;

% Calibration Maps here
omega = 0:0.5:200;  % Off-resonance,  B0
fam   = 1.0;    % Flip Angle Map, B1

for jj = 1:length(omega)

  % Generate Noise realizations
  for ii = 1:N
    % -- Generate Signal Curves Using CPU --
    [data_spgr(ii,:) data_ssfp_0(ii,:)  ] = mcdespot_sim_v2([1000 1 1/t1_f 1/t2_f 1./t1_m 1./t2_m fract tau omega(jj)*2*pi], alpha_spgr, alpha_ssfp, tr_spgr, tr_ssfp, SNR);
    [x               data_ssfp_180(ii,:)] = mcdespot_sim_v2([1000 1 1/t1_f 1/t2_f 1./t1_m 1./t2_m fract tau (omega(jj)*tr_ssfp*2*pi + pi)/tr_ssfp], alpha_spgr, alpha_ssfp, tr_spgr, tr_ssfp, SNR);
  end

  [fv rnrm] = cpMCDESPOT_model_fit(data_spgr, data_ssfp_0, data_ssfp_180, alpha_spgr, alpha_ssfp, tr_spgr, tr_ssfp, fam*ones([N 1]), zeros([N 1]), mode, debug);

  
  fv_omega_mean(jj,:) = mean(fv);
  fv_omega_std(jj,:)  = std(fv);
  fv = [];
end

toc;