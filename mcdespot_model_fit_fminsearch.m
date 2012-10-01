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
% v4.2 09-Feb-2011
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

function [fv rnrm] = mcdespot_model_fit(data_spgr, data_ssfp_0, data_ssfp_180, alpha_spgr, alpha_ssfp, tr_spgr, tr_ssfp, fam, omega, ig, debug)

tic;

%% O. Optimization Settings
% Nelder-Mead Downhill Simplex
optim = optimset('fminsearch');
optim.TolFun  = 1e-4;
optim.TolX    = 1e-4;
optim.Display = 'none';

% Apply fmap B1 correction to flip angles
alpha_spgr = fam * alpha_spgr;
alpha_ssfp = fam * alpha_ssfp;

% Preallocate outputs
fv   = zeros([size(data_spgr, 1)  6]);
rnrm = zeros([size(data_spgr, 1)  1]);

% Determine number of points
npts = size(find(~(sum(data_spgr, 2) == 0)),1);
pt   = 0;
fprintf('mcDESPOT Fitting:');

%% Loop over non-zero voxels in the image
for ii = find(~(sum(data_spgr, 2) == 0))'
  progressbar(pt/npts);
  pt = pt + 1;
  
  % Grab voxel data for SPGR & SSFP
  vox_data_spgr     = data_spgr(ii,:)'     ./ ig(ii,3);
  vox_data_ssfp_0   = data_ssfp_0(ii,:)'   ./ ig(ii,4);
  vox_data_ssfp_180 = data_ssfp_180(ii,:)' ./ ig(ii,4);
  
  % Grab flip angles for current voxel
  vox_alpha_spgr = alpha_spgr(ii,:)';
  vox_alpha_ssfp = alpha_ssfp(ii,:)';
  
  % Grab off-resonance value for this voxel
  vox_omega = omega(ii);
  
  % Initial guess for this voxel
  %         T1m          T1f          T2m          T2f           MWF   Tau
  vox_ig = [ig(ii,1)*0.8 ig(ii,1)*1.2 ig(ii,2)*0.8 ig(ii,3)*1.2  0.05  0.150];
  
  % Minimize with fminsearch
  [x res] = fminsearch(@mcdespot_model, vox_ig, optim);
  
  % Assign results to output vector
  fv(ii, 1:6) = x;
  rnrm(ii)    = res;

end % End voxels

disp('Done!');
toc;


%% Wrapper for C-Code
  function res = mcdespot_model(x)
    [resSPGR resSSFP_0 resSSFP_180] = cpMCDESPOT_residuals_SAH([x(1:6) vox_omega], vox_data_spgr, vox_data_ssfp_0, vox_data_ssfp_180, vox_alpha_spgr, vox_alpha_ssfp, tr_spgr, tr_ssfp, 1); % Extra 1 for tefix flag
    res = resSPGR + resSSFP_0 + resSSFP_180;
  end



end