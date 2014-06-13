% FUNCTION [] = run_despot1();
%
% FUNCTION to automate running DESPOT1-HIFI to generate flip angle map,
%          then T1 fitting routine to generate single-componant T1 values
%
% Inputs:
%    EXPLICIT: none
%    IMPLICIT: Start in same directory as _mcdespot_settings.mat
%
% Outputs:
%    EXPLICIT: None
%    IMPLICIT: New directory called singleComponent
%
% Samuel A. Hurley
% v5.1 05-Mar-2014
%
% Changelog:
%     v3.0 - Initial Version (using v3.0 to match other mcDESPOT commands)     (Sept-2010)
%     v3.1 - Test of SVN (Oct-2010)
%     v3.2 - Changed to 3D polynomial smoothing for input volumes (May-2011)
%     v3.3 - Normalizes input PD by max voxel signal in masked image. (Jul-2011)
%     v3.4 - Fixed normalization of 2nd (non-smoothed) T1 iteration   (Jan-2012)
%     v3.5 - Updated to be compatible with afi_flag & ideal_flag options (Feb-2012)
%     v3.6 - Fixed mask.spgr() bug when BET mask is not specified (Jan-2014)
%     v5.0 - Skipped v4 to make maj rev match between all mcDESPOT codes (Mar-2014)
%            Added support for Bloch-Siegert FA map
%     v5.1 - Use [dir.BASE dir.SPGR 'spgr_01'] instead of info_spgr for writing out
%            new NIfTI file headers. (Avoids FSL orientation mis-label if flags.reorient
%            is set. (Mar 2014)

function [] = run_despot1()

% Initiate Diary
diary('_mcdespot_log.txt');

% Display Banner
VER     = 5.1;
VERDATE = '05-Mar-2014';

% Display banner
disp('=== cpMCDESPOT - Multicomponent Relaxomtery Analysis ===');
disp('     DESPOT1-HIFI Script   (run_despot1)'                );
disp('     Samuel A. Hurley      shurley@wisc.edu'             );
disp(['    Version ' num2str(VER, '%01.1f') '         ' VERDATE]);
disp('     FOR USE ONLY AT UNIVERSITY OF WISCONSIN.'           );
disp('========================================================');

% Load in mcdespot settings file
load _mcdespot_settings;

% Check Flags
if ~isfield(flags, 'afi')
  afi_flag = 0;
else
  afi_flag = flags.afi;
end

% Starting time for DESPOT1
time.despot1_start = datetime();
disp(['DESPOT1 Started: ' time.despot1_start]);

% Define singleComponant location
dir.DESPOT1 = './singleComponent/';
mkdir(dir.DESPOT1);

% Load mcDESPOT Images
% load_mcdespot_series will choose the coreg data, if it exists
img = load_mcdespot_series();

% Free up memory
img.ssfp_0   = [];
img.ssfp_180 = [];

% Determine SPGR Normalization Factor
maxVal = max(img.spgr(:));
status.despot1_signalScale = maxVal;
disp(['Normalizing Data By Scale Factor: ' num2str(status.despot1_signalScale)]);

% Try to load the mask, if it exists
if isfield(status, 'mask') && status.mask == 1 
  % Load the mask
  mask = logical(load_nifti([dir.MASK status.maskname]));
  disp(['Using user-supplied mask.']);
else
  % Threshold above 0
  mask = img.spgr(:,:,:,1) > 0;
  disp('Using threshold mask.');
end

% Get size of data
dataSize = size(img.spgr);

%% AFI For FA Map
if afi_flag == 1
  %% Standard (non-VAFI) Mapping
  fam = afi_standard(img.afi_01, img.afi_02, 5) ./ alpha_afi;
  
  % 3D Smooth /w Mask
  disp('3D Smoothing AFI FA Map');
  fam_s = img_smooth_polyfit(fam,  mask, 3, 1);
  rnrm  = zeros(size(fam_s));

%% Bloch-Siegert for FA Map
elseif afi_flag == 2
  fam = img.afi_01 ./ (alpha_afi * 10);
  
  % 3D Smooth /w Mask
  disp('3D Smoothing AFI FA Map');
  fam_s = img_smooth_polyfit(fam,  mask, 3, 1);
  rnrm  = zeros(size(fam_s));
  
  
%% IR-SPGR (HIFI) for FA Map
else
  
  % Grab the necessary data
  spgr   = img.spgr;
  irspgr = img.irspgr;
  
  % Normalize data by maximum voxel value in masked SPGR image
  spgr   = spgr   ./ maxVal;
  irspgr = irspgr ./ maxVal;
  
%   % Smoothing & Masking raw data before fitting
%   disp('3D Smoothing data prior to fitting...');
%   diary('off');
%   spgr   = img_smooth_polyfit(spgr,   mask, 3, 1);
%   irspgr = img_smooth_polyfit(irspgr, mask, 3, 1);
%   diary('_mcdespot_log.txt');

  % DEBUG: Mask w/o Smoothing
  spgr = spgr .* repmat(mask, [1 1 1 length(alpha_spgr)]);
  
  %% Run DESPOT1-HIFI for B1 Map
  warning off %#ok<WNOFF>
  diary('off');
  opts.debug = 0;
  [pd r1 fam rnrm] = despot1hifi_model_fit(reshape(spgr, [dataSize(1)*dataSize(2)*dataSize(3) dataSize(4)]), alpha_spgr, tr_spgr, irspgr(:), alpha_irspgr, tr_irspgr, ti_irspgr/1000, npe_irspgr, opts); %#ok<ASGLU>
  warning on %#ok<WNON>
  diary('_mcdespot_log.txt');
  
  % Reshape Resulting FAM
  % Assume it is smooth enough from pre-smoothing fitted data
  fam_s = abs(reshape(fam, [dataSize(1) dataSize(2) dataSize(3)]));
  rnrm  =     reshape(rnrm, [dataSize(1) dataSize(2) dataSize(3)]);
end

%% Free up memory
img.afi_01 = [];
img.afi_02 = [];
img.irspgr = [];

%% Re-Grab and re-mask unsmoothed SPGR data
spgr   = img.spgr;
spgr   = spgr .* repmat(mask, [1 1 1 length(alpha_spgr)]);

% Normalize data by same value as smoothed SPGR
spgr   = spgr   ./ maxVal;

% Re-process using IRLS - 3 Iterations
warning off; %#ok<WNOFF>
diary('off');
[pd r1] = t1_fit_spgr_IRLS_afi(reshape(spgr, [dataSize(1)*dataSize(2)*dataSize(3) dataSize(4)]), alpha_spgr, tr_spgr, fam_s(:), 3);
warning on; %#ok<WNON>
diary('_mcdespot_log.txt');

% Rehape Results
r1  = reshape(r1,  [dataSize(1) dataSize(2) dataSize(3)]);
pd  = reshape(pd,  [dataSize(1) dataSize(2) dataSize(3)]);

% Save NIfTI
img_nifti_to_nifti(iminv(r1), [dir.BASE dir.SPGR 'spgr_01.nii'], [dir.DESPOT1 'DESPOT1-T1']);
img_nifti_to_nifti(fam_s,     [dir.BASE dir.SPGR 'spgr_01.nii'], [dir.DESPOT1 'DESPOT1-FAM']);
img_nifti_to_nifti(fam,       [dir.BASE dir.SPGR 'spgr_01.nii'], [dir.DESPOT1 'DESPOT1-FAM_Unsmooth']);
img_nifti_to_nifti(pd,        [dir.BASE dir.SPGR 'spgr_01.nii'], [dir.DESPOT1 'DESPOT1-PD']);
img_nifti_to_nifti(rnrm,      [dir.BASE dir.SPGR 'spgr_01.nii'], [dir.DESPOT1 'DESPOT1-Rnrm']);

% Plot the center slice
centerSlice = round(dataSize(3) / 2);

% T1 Plot
subplot(1,2,1);
imagesc(iminv(r1(:,:, centerSlice)), [0 3]);
colorbar;
colormap gray;
axis image;
axis off;
title 'T1 DESPOT1 [s]';

% FAM Plot
subplot(1,2,2);
imagesc(fam_s(:,:, centerSlice), [.5 1.5]);
colorbar;
colormap gray;
axis image;
axis off;
title 'FAM Pre-smoothed from HIFI';

% Save Figure
% savefig([dir.DESPOT1 'DESPOT1-Result.tif']);

% Done!
time.despot1_end = datetime();
disp(['DESPOT1-HIFI Complete: ' time.despot1_end]);

% DESPOT1 Flag
status.despot1 = 1;
  
save('_mcdespot_settings', 'tr_spgr', 'tr_irspgr', 'tr_ssfp', 'flags', ...
     'alpha_spgr', 'alpha_afi', 'alpha_irspgr', 'alpha_ssfp', 'ti_irspgr', 'npe_irspgr', ...
     'info_spgr', 'info_irspgr', 'info_afi', 'info_ssfp_0', 'info_ssfp_180', 'info_ideal', 'dir', 'time', 'status');
      
diary('off');

