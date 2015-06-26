%% Plot out 2D 'terrain' of mcDESPOT residuals
tic;

% Weighting for residuals
SPGRWEIGHT      = 2.75;
SSFPWEIGHT_HIGH = 2.15;
SSFPWEIGHT_LOW  = 0.45;

mwf  = 0.18;

DEBUG = 0;      % Debug mcDESPOT Fitting?

% Scan Parameters
tr_spgr = .0065;
tr_ssfp = .0050;

% Use flip angles from actual mcDESPOT v2 Protocol
alpha_spgr = [3 4 5 6 7 9 13 18];
alpha_ssfp = [12 16 21 27 33 40 51 68];

% ----   #1 - MWF Sweep, Hold Others Const  ----
t1_m = 0.465;
t1_f = 1.07;

t2_m = 0.026;
t2_f = 0.117;

tau  = 0.180;

% Calibration Maps here
omega   = 40.00;  % Off-resonance,  B0
fam     = 1.0;    % Flip Angle Map, B1
pd_spgr = 1;
pd_ssfp = 1;

% Parameter Vector
mcd_fv = [t1_m t1_f t2_m t2_f mwf tau omega pd_spgr pd_ssfp];

% Generate true signal
[T_spgr x x]              = sim_mcdespot_c(mcd_fv, alpha_spgr, tr_spgr, tr_ssfp, Inf); % No noise -- this is added elsewhere
[x T_ssfp_0 T_ssfp_180]   = sim_mcdespot_c(mcd_fv, alpha_ssfp, tr_spgr, tr_ssfp, Inf);

% Compute residuals over parameter ranges

% Ranges
t1_m_range = 0.300:.010:0.650;
t1_f_range = 0.500:.020:1.500;

t2_m_range = 0.001:0.0005:0.030;
t2_f_range = 0.050:.005:0.165;

mwf_range  = 0.01:0.005:0.30;
tau_range  = 0.025:0.010:0.600;

%% T1 Space
[X Y] = meshgrid(t1_m_range, t1_f_range);

sz = size(X(:));

x = [X(:) Y(:) repmat(t2_m, [sz 1]) repmat(t2_f, [sz 1]) repmat(mwf, [sz 1]) repmat(tau, [sz 1])];

% Compute residuals via CPU
resSPGR     = cpMCDESPOT_residuals_SAH(x', omega,  -1, T_spgr',     alpha_spgr', tr_spgr, 1);
resSSFP_0   = cpMCDESPOT_residuals_SAH(x', omega,   0, T_ssfp_0',   alpha_ssfp', tr_ssfp, 1);
resSSFP_180 = cpMCDESPOT_residuals_SAH(x', omega, 180, T_ssfp_180', alpha_ssfp', tr_ssfp, 1);

%     % DEBUG -- Equal To All
%     res = resSPGR + resSSFP_0 + resSSFP_180;

% Compute residual based on weighting
off_res_range = 1/2/tr_ssfp; % Range of 0->off_res_range of omega values

% |0%| -- 180 Only -- |33%| -- 180>0 -- |50%| -- 0>180 -- |66%| -- 0 Only -- |1/(2*omega)|

if omega < off_res_range * 0.33
  % Use SSFP-180 Only (SSFP-0 is zero signal)
  res = SPGRWEIGHT*resSPGR + (SSFPWEIGHT_HIGH+SSFPWEIGHT_LOW)*resSSFP_180;
elseif omega < off_res_range * 0.50
  % Use both, weight SSFP-180 Higher
  res = SPGRWEIGHT*resSPGR + SSFPWEIGHT_HIGH*resSSFP_180 + SSFPWEIGHT_LOW*resSSFP_0;
elseif omega < off_res_range * 0.66
  % Use both, weight SSFP-0 Higher
  res = SPGRWEIGHT*resSPGR + SSFPWEIGHT_HIGH*resSSFP_0   + SSFPWEIGHT_LOW*resSSFP_180;
else
  % Use only SSFP-0, since SS 0FP-180 is almost zero in this region
  res = SPGRWEIGHT*resSPGR + (SSFPWEIGHT_HIGH+SSFPWEIGHT_LOW)*resSSFP_0;
end

res_t1 = res;

res_t1 = reshape(res_t1, size(X));
surf(X, Y, res_t1);
xlabel 'T1_m [s]';
ylabel 'T1_f [s]';
zlabel 'Residual';
hold on;
plot(t1_m, t1_f, 'w+', 'MarkerSize',10, 'LineWidth', 2);
hold off;
zlim([0 .02]);
caxis([0 .02]);

%% T2 Space
[X Y] = meshgrid(t2_m_range, t2_f_range);

sz = size(X(:));

x = [repmat(t1_m, [sz 1]) repmat(t1_f, [sz 1]) X(:) Y(:) repmat(mwf, [sz 1]) repmat(tau, [sz 1])];

% Compute residuals via CPU
resSPGR     = cpMCDESPOT_residuals_SAH(x', omega,  -1, T_spgr',     alpha_spgr', tr_spgr, 1);
resSSFP_0   = cpMCDESPOT_residuals_SAH(x', omega,   0, T_ssfp_0',   alpha_ssfp', tr_ssfp, 1);
resSSFP_180 = cpMCDESPOT_residuals_SAH(x', omega, 180, T_ssfp_180', alpha_ssfp', tr_ssfp, 1);

%     % DEBUG -- Equal To All
%     res = resSPGR + resSSFP_0 + resSSFP_180;

% Compute residual based on weighting
off_res_range = 1/2/tr_ssfp; % Range of 0->off_res_range of omega values

% |0%| -- 180 Only -- |33%| -- 180>0 -- |50%| -- 0>180 -- |66%| -- 0 Only -- |1/(2*omega)|

if omega < off_res_range * 0.33
  % Use SSFP-180 Only (SSFP-0 is zero signal)
  res = SPGRWEIGHT*resSPGR + (SSFPWEIGHT_HIGH+SSFPWEIGHT_LOW)*resSSFP_180;
elseif omega < off_res_range * 0.50
  % Use both, weight SSFP-180 Higher
  res = SPGRWEIGHT*resSPGR + SSFPWEIGHT_HIGH*resSSFP_180 + SSFPWEIGHT_LOW*resSSFP_0;
elseif omega < off_res_range * 0.66
  % Use both, weight SSFP-0 Higher
  res = SPGRWEIGHT*resSPGR + SSFPWEIGHT_HIGH*resSSFP_0   + SSFPWEIGHT_LOW*resSSFP_180;
else
  % Use only SSFP-0, since SS 0FP-180 is almost zero in this region
  res = SPGRWEIGHT*resSPGR + (SSFPWEIGHT_HIGH+SSFPWEIGHT_LOW)*resSSFP_0;
end

res_t2 = res;

res_t2 = reshape(res_t2, size(X));
surf(X, Y, res_t2);
xlabel 'T2_m [s]';
ylabel 'T2_f [s]';
zlabel 'Residual';
hold on;
plot(t2_m, t2_f, 'w+', 'MarkerSize',10, 'LineWidth', 2);
hold off;
zlim([0 .02]);
caxis([0 .02]);


%% MWF/TAU Space
[X Y] = meshgrid(mwf_range, tau_range);

sz = size(X(:));

x = [repmat(t1_m, [sz 1]) repmat(t1_f, [sz 1]) repmat(t2_m, [sz 1]) repmat(t2_f, [sz 1])  X(:) Y(:)];

% Compute residuals via CPU
resSPGR     = cpMCDESPOT_residuals_SAH(x', omega,  -1, T_spgr',     alpha_spgr', tr_spgr, 1);
resSSFP_0   = cpMCDESPOT_residuals_SAH(x', omega,   0, T_ssfp_0',   alpha_ssfp', tr_ssfp, 1);
resSSFP_180 = cpMCDESPOT_residuals_SAH(x', omega, 180, T_ssfp_180', alpha_ssfp', tr_ssfp, 1);

%     % DEBUG -- Equal To All
%     res = resSPGR + resSSFP_0 + resSSFP_180;

% Compute residual based on weighting
off_res_range = 1/2/tr_ssfp; % Range of 0->off_res_range of omega values

% |0%| -- 180 Only -- |33%| -- 180>0 -- |50%| -- 0>180 -- |66%| -- 0 Only -- |1/(2*omega)|

if omega < off_res_range * 0.33
  % Use SSFP-180 Only (SSFP-0 is zero signal)
  res = SPGRWEIGHT*resSPGR + (SSFPWEIGHT_HIGH+SSFPWEIGHT_LOW)*resSSFP_180;
elseif omega < off_res_range * 0.50
  % Use both, weight SSFP-180 Higher
  res = SPGRWEIGHT*resSPGR + SSFPWEIGHT_HIGH*resSSFP_180 + SSFPWEIGHT_LOW*resSSFP_0;
elseif omega < off_res_range * 0.66
  % Use both, weight SSFP-0 Higher
  res = SPGRWEIGHT*resSPGR + SSFPWEIGHT_HIGH*resSSFP_0   + SSFPWEIGHT_LOW*resSSFP_180;
else
  % Use only SSFP-0, since SS 0FP-180 is almost zero in this region
  res = SPGRWEIGHT*resSPGR + (SSFPWEIGHT_HIGH+SSFPWEIGHT_LOW)*resSSFP_0;
end

res_mwf = res;

res_mwf = reshape(res_mwf, size(X));
surf(X, Y, res_mwf);
xlabel 'MWF';
ylabel 'Tau [s]';
zlabel 'Residual';
hold on;
plot(mwf, tau, 'w+', 'MarkerSize',10, 'LineWidth', 2);
hold off;
zlim([0 .004]);
caxis([0 .02]);
