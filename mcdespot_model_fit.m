% FUNCTION [fv rnrm] = mcdespot_model_fit(data_spgr, data_ssfp_0,
% data_ssfp_180, alpha_spgr, alpha_ssfp, tr_spgr, tr_ssfp, fam, omega, ig, debug)
%
% Fits SPGR & SSFP MRI signal to a two-component T1 & T2 model using the
% mcDESPOT technique.  Assumes acquisition @3T Scanner
%
%
% Inputs:
%          data_spgr     [NpxNf]  - SPGR data, voxels in columns, flip angles in rows
%          data_ssfp_0   [NpxNf]  - SSFP with 0 degree phase-cycling
%          data_ssfp_180 [NpxNf]  - SSFP with 180 degree phase-cycling
%          alpha_spgr, alpha_ssfp - flip angles, in degrees
%          tr_spgr, tr_ssfp       - TR of sequnes, in seconds
%          fam                    - relative error in B1
%          omega                  - off-resonance of B0 in Hz
%          ig                     - Np x 4 initial guess variable [t1 t2 pd_spgr pd_ssfp]
%
% Outputs:
%           fv [Np x 7] = [T1m T1f T2m T2f Fm Tau_m Omega]
%           rnrm        - residual norm
%
% Samuel A. Hurley
% University of Wisconsin
% v4.3 27-Jan-2012
%
% Chagelog:
%        v1.0 20-Oct-2009 - Initial Relase
%        v1.1 8-Dec-2009  - 
%        v2.0 30-Mar-2010 - Change weigihting scheme, initial guess values
%        v3.0 16-Jun-2010 - Add explicit fitting for proton density
%                           Changed function name to mcdespot_model_fit
%        v4.0 21-Oct-2010 - Trying SSFP model from Scheffler's Group
%        v4.1 01-Nov-2010 - Cleaned up zero voxel testing
%        V4.2 09-Feb-2011 - Code re-vamp based on new cpMCDESPOT model
%                           (fixed TE=TR/2 error)
%        V4.3 16-Oct-2011 - Now outputs 6 parameters (since PD and Omega's all work from singleComponant fitting now).
%                           Reverted from fminsearch back to Sean's algorithm.
%        V4.4 27-Jan-2012 - Added FMINSEARCH as 2nd step in fitting, using
%                           results from GaussianContraction as initial guess.

function [fv rnrm] = mcdespot_model_fit(data_spgr, data_ssfp_0, data_ssfp_180, alpha_spgr, alpha_ssfp, tr_spgr, tr_ssfp, fam, omega, ig, debug)

tic;

%% O. Define constants for weighting data, iterations, etc
RETAINED_TOP_RESIDUALS = 40;                       % Base number of solutions
TOP_SOLUTION           = 5;                        % Keep top 5 answers as solution

CONTRACTION_STEPS      = 4;
NUM_RANDOM_WALKS       = 28;
NUM_SAMPLES            = 128;  % Must be 128 to match #threads/kernel on GPU

% SIGNALSCALE = 1000;  % Hard-coded into GPU code
% MAX_ALPHA_SPGR = 40; % Hard-coded into GPU code 
% MAX_ALPHA_SSFP = 40; % Hard-coded into GPU code

% Weighting for residuals
SPGRWEIGHT      = 2.75;
SSFPWEIGHT_HIGH = 2.15;
SSFPWEIGHT_LOW  = 0.45;

% Nelder-Mead Downhill Simplex (fminsearch)
optim=optimset('fminsearch');
optim.TolFun = 1e-4;
optim.TolX   = 1e-4;
optim.MaxFunEvals = 100;
optim.Display = 'off';

% Counter for debug
debug_ctr = 0;

% Apply fmap B1 correction to flip angles
alpha_spgr = fam * alpha_spgr;
alpha_ssfp = fam * alpha_ssfp;

% Preallocate outputs
fv   = zeros([size(data_spgr, 1)  6]);
rnrm = zeros([size(data_spgr, 1)  1]);

% Determine number of points
npts = size(find(~(sum(data_spgr, 2) == 0)),1);
pt   = 0;

if debug == 0
  fprintf('mcDESPOT Fitting:');
else
  disp('mcDESPOT Fitting:');
end

%% Loop over non-zero voxels in the image
for ii = find(~(sum(data_spgr, 2) == 0))'
  
  if debug == 0
    % Only show progressbar when not doing debugging
    progressbar(pt/npts);
  end
  pt = pt + 1;
    
  % Grab SPGR & SSFP Data For the current voxel, normalize by proton density
  % *Transpose data because C-routine indexies via column vectors
  warning off;
  vox_data_spgr     = data_spgr(ii,:)'     ./ ig(ii,3);
  vox_data_ssfp_0   = data_ssfp_0(ii,:)'   ./ ig(ii,4);
  vox_data_ssfp_180 = data_ssfp_180(ii,:)' ./ ig(ii,4);
  warning on;

  % Grab flip angles for current voxel
  vox_alpha_spgr    = alpha_spgr(ii,:)';
  vox_alpha_ssfp    = alpha_ssfp(ii,:)';

  % Reset guess to initial guess
  guess = initialGuess3T();
  
  % Set the B0 to the DESPOT2-FM Omega
  vox_omega = omega(ii);

  %% Iterate over contraction steps
  for jj = 1:CONTRACTION_STEPS            % Loop for iterations
    
    % Generate NUM_RANDOM_WALKS * NUM_SAMPLES
    x = generatePoints(guess, NUM_RANDOM_WALKS*NUM_SAMPLES);
    
    % Compute residuals via CPU
    [resSPGR resSSFP_0 resSSFP_180] = cpMCDESPOT_residuals_SAH(x', vox_omega, vox_data_spgr, vox_data_ssfp_0, vox_data_ssfp_180, vox_alpha_spgr, vox_alpha_ssfp, tr_spgr, tr_ssfp, 1);
       
%     % DEBUG -- Equal To All 
%     res = resSPGR + resSSFP_0 + resSSFP_180;
    
    % Compute residual based on weighting
    off_res_range = 1/2/tr_ssfp; % Range of 0->off_res_range of omega values
    
    % |0%| -- 180 Only -- |33%| -- 180>0 -- |50%| -- 0>180 -- |66%| -- 0 Only -- |1/(2*omega)|
    
    if vox_omega < off_res_range * 0.33
      % Use SSFP-180 Only (SSFP-0 is zero signal)
      res = SPGRWEIGHT*resSPGR + (SSFPWEIGHT_HIGH+SSFPWEIGHT_LOW)*resSSFP_180;
    elseif vox_omega < off_res_range * 0.50
      % Use both, weight SSFP-180 Higher
      res = SPGRWEIGHT*resSPGR + SSFPWEIGHT_HIGH*resSSFP_180 + SSFPWEIGHT_LOW*resSSFP_0;
    elseif vox_omega < off_res_range * 0.66
      % Use both, weight SSFP-0 H 0ig 0her
      res = SPGRWEIGHT*resSPGR + SSFPWEIGHT_HIGH*resSSFP_0   + SSFPWEIGHT_LOW*resSSFP_180;
    else
      % Use only SSFP-0, since SS 0FP-180 is almost zero in this region
      res = SPGRWEIGHT*resSPGR + (SSFPWEIGHT_HIGH+SSFPWEIGHT_LOW)*resSSFP_0;
    end
    
%    % DEBUG: Show residual vectors
%     plot(1:(NUM_RANDOM_WALKS*NUM_SAMPLES), resSPGR, 1:(NUM_RANDOM_WALKS*NUM_SAMPLES), resSSFP_0, 1:(NUM_RANDOM_WALKS*NUM_SAMPLES), resSSFP_180);
%     ylim([0 .5]);
%     title(num2str(jj));
%     legend('SPGR', 'SSFP-0', 'SSFP-180');
%     pause(.5);
    
    % Sort residuals low->high
    [res idx] = sort(res);
    
    % Do not need to do if on last iteration
    if (jj < CONTRACTION_STEPS)
      % Update mean, min, max, stdev based on top residuals
      idx = idx(1:RETAINED_TOP_RESIDUALS);
      guess = updateGuess(x(idx,:));
    end
    
  end % End contraction steps
  
  % Return result as mean of top 5 solutions
  idx   = idx(1:TOP_SOLUTION);
  guess = updateGuess(x(idx,:));
  
  %% OPTIMIZATION STEP II: Use guess as initial guess for local optim (fminsearch)
  
%   % Gaussian Contraction + FMINSEARCH
%   [x rnrm(ii)]   = fminsearch(@mcdespot_model, guess(:,1)', optim);

  % DEBUG: Gaussian Contraction Only
  x = guess(:,1)';
  rnrm(ii) = res(1);
  
  %% Output vector is fv, also return mean of top 5 residuals
  fv(ii, :)   = x;
  
  %% DEBUG: Plot final fitted curve
  
  if debug == 1
    if mod(debug_ctr, 10) == 0
      disp('----------|---------|-----------|---------|-----|--------|-------|---------');
      disp('T1 Myelin | T1 Free | T2 Myelin | T2 Free | MWF |  Tau   |  Res  |Progress');
      disp('----------|---------|-----------|---------|-----|--------|-------|---------');
    end
    debug_ctr = debug_ctr + 1;
    disp(['  ' num2str(fv(ii,1)*1000, '%03.0f') ' ms  | ' num2str(fv(ii,2)*1000, '%04.0f') ' ms |   ' num2str(fv(ii,3)*1000, '%02.0f') ' ms   |  ' num2str(fv(ii,4)*1000, '%03.0f') ' ms | ' num2str(fv(ii,5)*100, '%02.0f') '% | ' num2str(fv(ii,6)*1000, '%03.0f') ' ms | ' num2str(rnrm(ii), '%03.3f') ' | ' num2str(pt*100/npts) '%']);
  end
  
  if debug == 1
    % Setup Figure
    if ~exist('dbfig', 'var')
      dbfig = figure;
    end

    % Generate simulated curve based on fitting
    a = 1:70;
    [s_spgr s_ssfp_0 s_ssfp_180] = sim_mcdespot(fv(ii,1:6), vox_omega, a, tr_spgr, tr_ssfp, Inf);

    % Plot actual vs fitted data
    figure(dbfig);
    plot(vox_alpha_spgr, vox_data_spgr, 'ok', vox_alpha_ssfp, vox_data_ssfp_0, 'or', vox_alpha_ssfp, vox_data_ssfp_180, 'ob', a, s_spgr, '-k', a, s_ssfp_0, '-r', a, s_ssfp_180, '-b');

    % Figure captions & Legend
    title(['mcDESPOT Signal Curves, RES=' num2str(rnrm(ii))]);
    % legend('MRI SPGR', 'MRI SSFP-0', 'MRI SSFP-180', 'Fit SPGR', 'Fit SSFP-0', 'Fit SSFP-180');
    xlabel('Flip Angle [degrees]');
    ylabel('Signal [a.u.]');
    drawnow;
    
  end

end % End voxels

% Close progressbar
progressbar(1);
toc;

%% ------------------------ Helper Functions Below -------------------------

% Wrapper for fminsearch call of mcDESPOT model
  function res = mcdespot_model(x)
    [resSPGR resSSFP_0 resSSFP_180] = cpMCDESPOT_residuals_SAH([x(:,1), x(:,2), x(:,3), x(:,4), x(:,5), x(:,6)]',vox_omega(:), vox_data_spgr, vox_data_ssfp_0, vox_data_ssfp_180, vox_alpha_spgr, vox_alpha_ssfp, tr_spgr, tr_ssfp, 1); %#ok<ASGLU>
    
    % DEBUG: Completly Equal Wt
    res = resSPGR + resSSFP_0 + resSSFP_180;
    
%     if vox_omega < off_res_range * 0.33
%       % Use SSFP-180 Only (SSFP-0 is zero signal)
%       res = SPGRWEIGHT*resSPGR + (SSFPWEIGHT_HIGH+SSFPWEIGHT_LOW)*resSSFP_180;
%     elseif vox_omega < off_res_range * 0.50
%       % Use both, weight SSFP-180 Higher
%       res = SPGRWEIGHT*resSPGR + SSFPWEIGHT_HIGH*resSSFP_180 + SSFPWEIGHT_LOW*resSSFP_0;
%     elseif vox_omega < off_res_range * 0.66
%       % Use both, weight SSFP-0 H 0ig 0her
%       res = SPGRWEIGHT*resSPGR + SSFPWEIGHT_HIGH*resSSFP_0   + SSFPWEIGHT_LOW*resSSFP_180;
%     else
%       % Use only SSFP-0, since SS 0FP-180 is almost zero in this region
%       res = SPGRWEIGHT*resSPGR + (SSFPWEIGHT_HIGH+SSFPWEIGHT_LOW)*resSSFP_0;
%     end
  end

%% Initial Guess for 3.0T Magnet, human brain
  function ig = initialGuess3T()
    %         [mean  min  max  std]
    % T1m
    ig(1,:) = [0.50 0.30 0.70 0.05];
    % T1f
    ig(2,:) = [1.40 0.80 2.80 0.200];

    % T2m
    ig(3,:) = [.015 .001 .040 .005];
    % T2f
    ig(4,:) = [.080 .050 .140 .020];

    % MWF
    ig(5,:) = [.150 .001 .250 .100];
    % Tau
    ig(6,:) = [.120 .075 .250 .025];
  end

%% Generate guess points based on gaussian model, from initial guess

% % DEBUG DEBUG -- make all pts the same
% % Based on Gaussian distributions
%   function pts = generatePoints(ig, nguess)
% 
%     % Preallocate x
%     pts = zeros([nguess 6]);
% 
%     % For each paramter
%     for  kk = 1:6
%       % Mean + (random gaussian)*std
%       tmp = ones([nguess 1])*ig(kk,1);
% 
%       pts(:,kk) = tmp;
%     end
% 
%   end

% Based on Gaussian distributions
  function pts = generatePoints(ig, nguess)

    % Preallocate x
    pts = zeros([nguess 6]);

    % For each paramter
    for  kk = 1:6
      % Mean + (random gaussian)*std
      tmp = ones([nguess 1])*ig(kk,1) + randn([nguess 1])*ig(kk,4);
      % Set all values below minimum to minimum
      tmp((tmp < ig(kk,2))) = ig(kk,2);
      % Set all values above maximum to maximum
      tmp((tmp > ig(kk,3))) = ig(kk,3);

      pts(:,kk) = tmp;
    end

  end

%% Update guess by calculating mean, min, max, and std of new parameters
  function g = updateGuess(xvals)

    g = zeros([6 4]);

    % Loop over each parameter
    for kk = 1:6
      % Update guess matrix with calcuated mean, min, max, and std
      g(kk,:) = [mean(xvals(:,kk)) 0.95*min(xvals(:,kk)) 1.05*max(xvals(:,kk)) std(xvals(:,kk))];
    end

  end

end
