% Function [pd r1 fam res] = despot1hifi_model_fit(data, alpha, tr, data_ir, alpha_ir, tr_ir, ti_ir, npe, opts)
%
% Function to fit linearized SPGR signal to VFA/DESPOT1-HIFI model
%
% Inputs:
%   data  - (Num Points) x (Num of Flip Angles) SPGR Data
%   alpha - a vector of flip angles [degrees]
%   tr    - repitition time [s]
%
%   data_ir - (Num Points) x (Num of Flip Angles/TI) IR-SPGR Data
%   alpha_ir - a vector of IR-SPGR flip angles [degrees]
%   tr_ir    - IR-SPGR repitition time [s]
%   ti_ir    - IR-SPGR inversion time  [s]
%   npe      - Number of phase encodes between inversion pulses.
%
%   opts  - a structure containing optional settings
%       points   - use a subset of points for 
%       debug    - debug by plotting data pts and fit curve
% 
%
% Outputs:
%   pd    - Proton Density [a.u.]
%   r1    - Spin-lattice relaxation rate [s^-1]
%   fam   - B1 Error (Relative to 1)
%   res   - Sum-of-squares residual
%
% References:
%   DESPOT1:     Deoni, SCL et al MRM 2003.
%
%   DESPOT-HIFI: Deoni, SCL. JMRI 2007;26:1106-1111
%
% Samuel A. Hurley
% University of Wisconsin
% v1.0 14-Sept-2010
%
% Chagelog:
%      v1.0 - Initial version based on t1_fit_spgr_IRLS (A.A. Samsonov, May 2005)
%                                  and t1_fit_ir_signa  (S.A. Hurley,   May 2010)

function [pd r1 fam rnrm] = despot1hifi_model_fit(data, alpha, tr, data_ir, alpha_ir, tr_ir, ti_ir, npe, opts)

tic;

% OO. Optomization settings
% Nelder-Mead Downhill Simplex (fminsearch)
optim=optimset('fminsearch');
optim.TolFun = 1e-4;
optim.TolX   = 1e-4;
optim.Display = 'off';

% Trust Region-Reflective (lsqnonlin)
% optim=optimset('lsqnonlin');
% optim.TolX   = 1e-4;
% optim.Display = 'off';

% O. Constants
% Seans 3.0T PSD -- MODE 1
% TI_SCALE       = 0.84;
% PD_SCALE       = 0.975;   % SD - Scale of PD for IR-SPGR rel to SPGR
% READOUT_PULSES = npe + 2;
% EFF            = -1 + cos(0.97 * pi);    % Inversion effiency (fixed)

% Seans 3.0T PSD -- MODE2
TI_SCALE       = 0.90;
PD_SCALE       = 0.975;                  % SD - Scale of PD for IR-SPGR rel to SPGR
READOUT_PULSES = npe/2 + 2;
EFF            = -1 + cos(0.97 * pi);    % Inversion effiency (fixed)


% I.Check input arguments
switch nargin
  case 8
    % No options struct specified, use default options
    opts   = struct();    % Use default options
  case 9
    if ~isstruct(opts)
      error('Opts must be a structure.  Type ''help despot1hifi_model_fit'' for more information');
    end
  otherwise
    error('You must supply 7 or 8 input arguments.  Type ''help despot1hifi_model_fit'' for more information');
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

% Banner
disp('==============DESPOT1-HIFI T1===================')

% II. Check data dimensions
if size(alpha, 2) ~= size(data, 2)
  error('Number of SPGR flip angles and supplied SPGR data points do not match.');
end

if size(alpha_ir, 2) ~= size(data_ir, 2)
  error('Number of IR-SPGR flip angles and supplied IR-SPGR data points to not match');
end

% Determine which points to use
if ~exist('points', 'var')
  % If not specified, use all points
  points = 1:length(alpha);
else
  % Use specified points
  if length(points) > length(alpha)
    error('You specifed more fitting points than avalible flip angles');
  end
  
  alpha = alpha(points);
  disp('NOTE: Subset of SPGR flip angles used for fitting.');
  
  % Grab only the flip angles specified
  data = data(:,points);
end

disp(['SPGR Flip Angles: ' num2str(alpha, '%01.0f ')]);
disp(['IR   Flip Angles: ' num2str(alpha_ir, '%01.0f ')]);

% Convert to radians
alpha    = alpha .* pi/180;
alpha_ir = alpha_ir .* pi/180;

% Preallocate output [pd r1 fam]
pd   = zeros([size(data,1) 1]);
r1   = zeros([size(data,1) 1]);
fam  = zeros([size(data,1) 1]);
rnrm = zeros([size(data,1) 1]);


% Number of non-zero data points
npts = size(find(~(sum(data, 2) == 0)),1);
pt   = 0;
fprintf('DESPOT1-HIFI Fitting:');
if debug > 0
  % Padding for twaitbar
  fprintf('                             ');
end

% III. Perform DESPOT1-HIFI Fit in Each Voxel
for ii = find(~(sum(data, 2) == 0))';
  % Grab voxel data
  vox_data    = data(ii,:);
  vox_data_ir = data_ir(ii,:);
  
  % Update Progress
  if debug == 0
    progressbar(pt/npts);
  else
    twaitbar(pt/npts);
  end
  pt = pt + 1;
  
  % Initial Guess - B1 correction of 1
  x = [1]; %#ok<NBRAK>
  
  % Minimize /w Levenburg-Marquardt
  [x res] = fminsearch(@despot1hifi_model, x, optim);
  
  % Save results
  pd(ii)   = pd_mod;
  r1(ii)   = r1_mod;
  fam(ii)  = x(1);
  rnrm(ii) = res;
  
  % PLOT FOR DEBUGGING FINAL FIT RESULT
  if debug == 1 && (toc - dbgtime) > 1
    % Only re-draw every 1 second
    dbgtime = toc;
    
    % Setup Figure
    if ~exist('dbfig', 'var')
      dbfig = figure;
    end
    
    % Generate Sim Data
    old_alpha = alpha;
    sim_alpha = (1:1:40) .* pi/180;
    alpha = sim_alpha;
    sim_data = despot1hifi_model(x,1);        % The 1 flags debug mode in the residual
    alpha = old_alpha;                        % Leaves alpha vector untouched
    
    % Plot actual vs fitted data
    figure(dbfig);
    plot(alpha*180/pi, data(ii,:), 'ob', sim_alpha*180/pi, sim_data(1:(end-1)), '-b', alpha_ir*180/pi, data_ir(ii,:)*4, 'or', alpha_ir*180/pi, sim_data(end)*4, 'xr');
    
    % Figure captions & Legend
    title(['B1 Error: ' num2str(x, '%0.3f') ' PD: ' num2str(pd_mod, '%0.0f') ' T1: ' num2str(1/r1_mod, '%0.2f') ' s']);
    legend('MRI SPGR', 'Fit SPGR', 'MRI IR-SPGR', 'Fit IR-SPGR');
    xlabel('Flip Angle [degrees]');
    ylabel('Signal [a.u.]');
    drawnow;
  end
  
end

% Done.
progressbar(1);
toc;

% IV. DESPOT1-HIFI Model Function
  function res = despot1hifi_model(fam_guess,  dbgmode) %#ok<INUSD>
    
    % Flip angle terms
    sina = sin(alpha.*fam_guess);
    cosa = cos(alpha.*fam_guess);
    tana = tan(alpha.*fam_guess);
    
    % Compute PD and R1 with LLS
    if ~exist('dbgmode', 'var')
      X = vox_data ./ tana;
      Y = vox_data ./ sina;
      
      p         = polyfit_fast(X,Y,1);
      slope     = p(1);
      intercept = p(2);
      
      pd_mod = intercept/(1-slope);
      r1_mod = -log(slope)/tr;
    end
    
    % SPGR Residual
    E1   = exp(-tr.*r1_mod);
    spgr = pd_mod.*(1-E1).*sina ./ (1-E1.*cosa);
    
    % IR-SPGR Residual
    ti         = ti_ir .* TI_SCALE;      % Re-scale TI time, based on Sean's Code
    full_tr_ir = ti + READOUT_PULSES*tr_ir;
    irspgr     = PD_SCALE .* pd_mod .* sin(alpha_ir.*fam_guess) .* (1 + EFF .* exp(-ti.*r1_mod) + exp(-full_tr_ir.*r1_mod));
    irspgr     = abs(irspgr);
    
    % Output Residual
    if ~exist('dbgmode', 'var')
      res = norm(spgr-vox_data) + norm(irspgr-vox_data_ir);  % Sum-of-squares residual
    else
      res = [spgr irspgr];
    end
    
  end


end
