% Function [pd r2 omega rnrm] = despot2fm_model_fit(data_0, data_180, alpha, tr, r1, pd, fam, opts)
%
% Function to fit single-componant DESPOT2 model, including off-resonance effects
% and extra term to account for TE
%
% NOTE: Only requires data_0 to be masked
%
% Inputs:
%   data_0    - (Num Points) x (Num of Flip Angles) SSFP /w   0 phase cycling
%   data_180  - (Num Points) x (Num of Flip Angles) SSFP /w 180 phase cycling
%   alpha     - a vector of flip angles [degrees]
%   tr        - repitition time [s]
%   r1        - Spin-lattice relaxation map, as fit from DESPOT1
%   pd        - Proton density relaxation map, as fit from DESPOT1, for initial guess
%   fam       - flip angle map to account for B1 inhomogeniety
%
%   opts  - a structure containing optional settings
%       points   - use a subset of points for 
%       debug    - 0 = debug by plotting data pts and fit curve
%                  1 = debug by showing GA fitness plots & plot of curve
%
% Outputs:
%   pd    - Proton density [a.u.]
%   r2    - Spin-spin relaxation rate [s^-1]
%   omega - B0 off-resonance [Hz]
%   rnrm  - Sum-of-squares residual
%
% References:
%   DESPOT2-FM: 
%
% Samuel A. Hurley
% University of Wisconsin
% v5.0 14-Nov-2012
%
% Chagelog:
%      v1.0 - Initial version based on despot1hifi_model_fit     (Feb-2011)
%      v2.0 - Switched to GA (Genetic Algorithm) fitting routine (Jun-2011)
%      v2.1 - Normalize by mean of MR signal to avoid explicit fitting for PD in GA, similar to Sean's mcDESPOT (Jun-2011)
%      v2.2 - Two-pass process to fit first for Omega ig (assuming fixed T2), 
%             then again for T2, PD, and Omega. Option for external B0 map  (Oct-2011)
%      v2.3 - Input PD initial guess from DESPOT1. Use R1 to formulate T2
%             initial guess based on Vasily constraints. (Nov-2011)
%      v2.4 - Removed unused GA stuff, updated pd initial guess for ext. supplied B0 map (Jan-2012)
%      v5.0 - Update to be compatible /w new cpDESPOT2 C-Code (Nov-2012)

function [pd r2 omega rnrm] = despot2fm_model_fit(data_0, data_180, alpha, tr, r1, pd_despot1, fam, opts)

tic;

% Constant: Assumption of ratio of T2 to T1 for 1st pass
T2_T1_RATIO = 0.066; % SHP brain ratio, Yarhykh assumption is 22.2, or 1/0.045

% Nelder-Mead Downhill Simplex (fminsearch)
optim=optimset('fminsearch');
optim.TolFun = 1e-4;
optim.TolX   = 1e-4;
optim.MaxFunEvals = 100;
optim.Display = 'off';

% I.Check input arguments
switch nargin
  case 7
    % No options struct specified, use default options
    opts   = struct();    % Use default options
  case 8
    if ~isstruct(opts)
      error('Opts must be a structure.  Type ''help despot2fm_model_fit'' for more information');
    end
  otherwise
    error('You must supply 7 or 8 input arguments.  Type ''help despot2fm_model_fit'' for more information');
end

% Grab Options or Set Defaults
if isfield(opts, 'points')
  points = opts.points;
end

if isfield(opts, 'debug')
  debug = opts.debug;
  if debug > 0
    dbgtime = toc;
  end
else
  % Default is to turn off debugging
  debug = 0;
end

if isfield(opts, 'ext_omega')
  % Externally supplied B0 map, in Hz.
  ext_omega = opts.ext_omega;
end


% Banner
disp('==============DESPOT2-FM T2===================');

% II. Check data dimensions
if size(alpha, 2) ~= size(data_0, 2)
  error('Number of SSFP-0 flip angles and supplied SPGR data points do not match.');
end

if size(alpha, 2) ~= size(data_180, 2)
  error('Number of SSFP-180 flip angles and supplied SPGR data points do not match.');
end

% Determine which points to use
if exist('points', 'var')
  % Use specified points
  if length(points) > length(alpha)
    error('You specifed more fitting points than avalible flip angles');
  end
  
  alpha = alpha(points);
  disp('NOTE: Subset of SSFP flip angles used for fitting.');
  
  % Grab only the flip angles specified
  data_0   = data_0(:,points);
  data_180 = data_180(:,points);
end

disp(['SSFP Flip Angles: ' num2str(alpha, '%01.0f ')]);

% Apply Flip Angle Correction
alpha = fam * alpha;

% Convert to radians
alpha    = alpha .* pi/180;

% Preallocate output [pd r2 omega]
pd    = zeros([size(data_0, 1) 1]);
r2    = zeros([size(data_0, 1) 1]);
omega = zeros([size(data_0, 1) 1]);
rnrm  = zeros([size(data_0, 1) 1]);

% Number of non-zero data points
npts = size(find(~(sum(data_0, 2) == 0)),1);
pt   = 0;
fprintf('DESPOT2-FM Fitting:');
if debug > 0
  % Padding for twaitbar
  fprintf('                             ');
end

% III. Perform DESPOT2-FM Fit in Each Voxel
for ii = find(~(sum(data_0, 2) == 0))';
  % Grab voxel data
  vox_data_0   = data_0(ii,:);
  vox_data_180 = data_180(ii,:);
  
  vox_alpha    = alpha(ii, :);
  vox_r1       = r1(ii);
  vox_pd_d1    = pd_despot1(ii);
  
  % Update Progress
  if debug == 0
    progressbar(pt/npts);
  else
    twaitbar(pt/npts);
  end
  pt = pt + 1;
 
  % Minimize with Fminsearch
  if ~exist('ext_omega', 'var');
    % Initial pass: fix T2, fit for initial guess of Omega
    % Fit for PD and Off-Res
    fix_flag = 1;
    ig = [vox_pd_d1 10];                   % PD / Omega initial guess
    r2_ig = vox_r1 / T2_T1_RATIO;          % Fix T2 based on T2/T1 ~= 0.045 (0.055 from Yarnykh NI2004)
    [x1] = fminsearch(@despot2fm_model, ig, optim);
  
    % Second pass - set omega IG from previous step
    fix_flag = 0;
    ig = [x1(1) vox_r1/T2_T1_RATIO x1(2)]; % PD from X1, R2 = R1/0.045, Omega from X1
    [x residual]   = fminsearch(@despot2fm_model, ig, optim);
    
  else
    % Single pass with fixed Omega (from ext B0 map)
    vox_omega = ext_omega(ii);             % Grab ith voxel of external B0 map.
    fix_flag  = 2;
    
    ig = [vox_pd_d1 vox_r1/T2_T1_RATIO];   % PD / T2 initial guess
    [x residual]   = fminsearch(@despot2fm_model, ig, optim);
    x(3) = vox_omega;
  end
  
  % Save results
  pd(ii)    = x(1);
  r2(ii)    = x(2);
  omega(ii) = x(3);
  rnrm(ii)  = residual;
  
  % PLOT FOR DEBUGGING FINAL FIT RESULT
  if debug > 0 && (toc - dbgtime) > .250
    % Only re-draw every 1 second
    dbgtime = toc;
    
    % Setup Figure
    if ~exist('dbfig', 'var')
      dbfig = figure;
    end
    
    % Generate Sim Data
    old_alpha = vox_alpha;
    sim_alpha = (1:1:90) .* pi/180;
    
    vox_alpha = sim_alpha;
    sim_data  = despot2fm_dbg(x,1);
    vox_alpha = old_alpha;                        % Leaves alpha vector untouched
    
    % Plot actual vs fitted data
    figure(dbfig);
    
    plot(vox_alpha*180/pi, data_180(ii,:), 'ob', vox_alpha*180/pi, data_0(ii,:), 'or', sim_alpha*180/pi, sim_data(1:(end/2)), '-b', sim_alpha*180/pi, sim_data((end/2 + 1):end), '-r');
    
    % Figure captions & Legend
    title([' PD: ' num2str(x(1), '%01.3f') ' T2: ' num2str(1/x(2), '%0.2f') ' s Omega: ' num2str(x(3))]);
    legend('SSFP-180 Data', 'SSFP-0 Data', 'Fit Curve 180', 'Fit Curve 0');
    xlabel('Flip Angle [degrees]');
    ylabel('Signal [a.u.]');
    drawnow;
    pause(.5);
  end
  
end

% Done.
progressbar(1);
toc;

% IV. DESPOT2-FM C-Function Wrapper
  function res = despot2fm_model(x)
    
    % Setup FV Call
    if fix_flag == 1
      fv = [x(1) 1./vox_r1 1./r2_ig  x(2)     ];
    elseif fix_flag == 0
      fv = [x(1) 1./vox_r1 1./x(:,2) x(:,3)   ];
    elseif fix_flag == 2
      fv = [x(1) 1./vox_r1 1./x(:,2) vox_omega];
    end
    
    % Call C-Code
    res =       cpDESPOT2_residuals_SAH(fv, 0,   vox_data_0',   vox_alpha', tr, 1)';
    res = res + cpDESPOT2_residuals_SAH(fv, 180, vox_data_180', vox_alpha', tr, 1)';

  end

% IV. DESPOT2-FM DEBUG M-Function
  function res = despot2fm_dbg(x, dbgmode) %#ok<INUSD>
    
    % Pull out model parameters
    pd_mod    = x(1);
    r2_mod    = x(2);
    omega_mod = x(3);
    
    % Flip angle terms
    sina  = sin(vox_alpha);
    cosa  = cos(vox_alpha);

    E1   = exp(-tr.*vox_r1);
    E2   = exp(-tr.*r2_mod);
    
    % SSFP-180 Signal
    beta = omega_mod*2*pi*tr + pi;
    sinb = sin(beta);
    cosb = cos(beta);
    
    denom = (1-E1.*cosa) .* (1-E2.*cosb) - E2.*(E1-cosa).*(E2-cosb);
    
    Mx = pd_mod .* ((1-E1).*E2.*sina.*sinb)      ./ denom;
    My = pd_mod .* ((1-E1).*E2.*sina.*(cosb-E2)) ./ denom;
    
    ssfp_180 = sqrt(Mx.^2 + My.^2) * sqrt(exp(tr.*r2_mod));
    
    % SSFP-0 Signal
    beta = omega_mod*2*pi*tr;
    sinb = sin(beta);
    cosb = cos(beta);
    
    denom = (1-E1.*cosa) .* (1-E2.*cosb) - E2.*(E1-cosa).*(E2-cosb);
    
    Mx = pd_mod .* ((1-E1).*E2.*sina.*sinb)      ./ denom;
    My = pd_mod .* ((1-E1).*E2.*sina.*(cosb-E2)) ./ denom;
    
    ssfp_0   = sqrt(Mx.^2 + My.^2) * sqrt(exp(tr.*r2_mod));
    
    % Output Residual
    if ~exist('dbgmode', 'var')
      res = norm([ssfp_180 ssfp_0] - [vox_data_180 vox_data_0]); % SOS Res
    else
      res = [ssfp_180 ssfp_0];
    end
    
  end

end

