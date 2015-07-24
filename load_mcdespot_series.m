% FUNCTION [img alpha_spgr alpha_ssfp] = load_mcdespot_series(mcdespot_settings)
%
% FUNCTION to automate loading of mcDESPOT raw images into Matlab variables
%
% Inputs:
%    Explicit: None
%    Implicit: Start in same directory as _mcdespot_settings.mat
%
% Outputs:
%    img  STRUCT of 4-D mcDESPOT images
%      img.spgr img.irspgr img.ssfp_0 img.ssfp_180
%
% Samuel A. Hurley
% Pouria Mossahebi
% University of Wisconsin
% University of Oxford
% v5.0 7-Jul-2015
%
% Changelog:
%     v1.0 - Based off of run_mcdespot (Mar-2010)
%     v3.0 - Changes to use the new ssfp_0 naming conventions (May-2011)
%     v3.1 - Eliminated stuff. (May-2010)
%     v3.2 - Supports datasets w/o IR-SPGR (ex: AFI acquisitions or Vnmr)  (Dec 2011)
%     v3.3 - Updated to be compatible with afi_flag & ideal_flag options   (Feb-2012)
%     v5.0 - Remove 'toc' at end. Update version number to match release   (Jul-2015)

function [img alpha_spgr alpha_ssfp] = load_mcdespot_series(mcdespot_settings) %#ok<STOUT>
tic;

% Load in mcdespot settings file
if ~exist('mcdespot_settings', 'var')
  load('_mcdespot_settings');
else
  load(mcdespot_settings);
end

% Check for flags
if ~isfield(flags, 'afi')
  afi_flag = 0;
else
  afi_flag = flags.afi;
end

if ~isfield(flags, 'ideal')
  ideal_flag = 0;
else
  ideal_flag = flags.ideal;
end

if isfield(status, 'coreg') && status.coreg == 1
  basedir = dir.COREG;
  disp('Loading Coregistered Data.');
else
  basedir = dir.BASE;
  disp('Loading Original Data.');
end

img_suf = '.nii';

% Only if IR-SPGR Dataset Exists (for AFI/Vnmr Datasets)
if afi_flag == 0
  fprintf('Loading IR-SPGR Data: ');
  progressbar(0);
  img.irspgr = load_nifti([basedir dir.IRSPGR   'irspgr'   img_suf]);
  progressbar(1);
  
  % Preallocate other data
  img.spgr     = zeros([size(img.irspgr) length(alpha_spgr)]);
  img.ssfp_0   = zeros([size(img.irspgr) length(alpha_ssfp)]);
  img.ssfp_180 = zeros([size(img.irspgr) length(alpha_ssfp)]);
else
  fprintf('Loading AFI Data: ');
  progressbar(0);
  img.afi_01 = load_nifti([basedir dir.EXT_CAL 'afi_01' img_suf]);
  progressbar(.5);
  img.afi_02 = load_nifti([basedir dir.EXT_CAL 'afi_02' img_suf]);
  progressbar(1);
end

% Load in data
fprintf('Loading SPGR Data:    ');
for ii = 1:length(alpha_spgr)
  progressbar(ii/length(alpha_spgr));
  img.spgr(:,:,:,ii)   = load_nifti([basedir dir.SPGR 'spgr_' num2str(ii, '%02.0f') img_suf]);
end

fprintf('Loading SSFP-0 Data:  ');
for ii = 1:length(alpha_ssfp)
  progressbar(ii/length(alpha_ssfp));
  img.ssfp_0(:,:,:,ii) = load_nifti([basedir dir.SSFP_0   'ssfp_0_' num2str(ii, '%02.0f') img_suf]);
end

fprintf('Loading SSFP-180 Data:');
for ii = 1:length(alpha_ssfp)
  progressbar(ii/length(alpha_ssfp));
  img.ssfp_180(:,:,:,ii) = load_nifti([basedir dir.SSFP_180 'ssfp_180_' num2str(ii, '%02.0f') img_suf]);
end

% Load in Omega
if ideal_flag == 1
  fprintf('Loading IDEAL B0 Data:');
  progressbar(0);
  img.omega = load_nifti([basedir dir.EXT_CAL 'IDEAL-Omega' img_suf]);
  progressbar(1);
end

disp('Done!');
