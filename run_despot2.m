% FUNCTION [] = run_despot2();
%
% FUNCTION to automate running DESPOT2-FM to generate simultaneous
%          single-componant T2 and B0 maps.
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
% v5.1 30-Apr-2014
%
% Changelog:
%     v3.0 - Initial Version (using v3.0 to match other mcDESPOT commands) (Jun-2011)
%     v3.1 - Updated to be compatible with afi_flag & ideal_flag options   (Feb-2012)
%     v5.0 - Skipped v4 to make maj rev match between all mcDESPOT codes (Mar-2014)
%     v5.1 - Use [dir.BASE dir.SPGR 'spgr_01'] instead of info_spgr for writing out
%            new NIfTI file headers. (Avoids FSL orientation mis-label if flags.reorient
%            is set.) (Apr-2014)



function [] = run_despot2()

% Initiate Diary
diary('_mcdespot_log.txt');

% Display Banner
VER     = 3.1;
VERDATE = '29-Feb-2012';

% Display banner
disp(['=== cpMCDESPOT - Multicomponent Relaxomtery Analysis ===']);
disp(['     DESPOT2-FM Script     (run_despot2)'                ]);
disp(['     Samuel A. Hurley      shurley@wisc.edu'             ]);
disp(['     Version ' num2str(VER, '%01.1f') '          ' VERDATE]);
disp(['     FOR USE ONLY AT UNIVERSITY OF WISCONSIN.'           ]);
disp(['========================================================']);

% Load in mcdespot settings file
load _mcdespot_settings;

% Check Flags
if ~isfield(flags, 'ideal')
  ideal_flag = 0;
else
  ideal_flag = flags.ideal;
end

% Starting time for coreg data
time.despot2_start = datetime();
disp(['DESPOT2-FM Started: ' time.despot2_start]);

% Define singleComponant location
dir.DESPOT1 = './singleComponent/';
mkdir(dir.DESPOT1);

% Load mcDESPOT Images
% load_mcdespot_series will choose the coreg data, if it exists
img = load_mcdespot_series();

% Grab the necessary data
ssfp_0   = img.ssfp_0;
ssfp_180 = img.ssfp_180;

% Free up memory
img.spgr   = [];
img.irspgr = [];
img.afi_01 = [];
img.afi_02 = [];

% Try to load single componant (T1 & FAM) maps as input to DESPOT2
if isfield(status, 'despot1') && status.despot1 == 1 %#ok<NODEF>
  % Load PD, T1, & FA Maps
  r1  = iminv(load_nifti([dir.DESPOT1 'DESPOT1-T1.nii' ]));
  pd  =       load_nifti([dir.DESPOT1 'DESPOT1-PD.nii' ]);
  fam =       load_nifti([dir.DESPOT1 'DESPOT1-FAM.nii']);
else
  % Cannot do DESPOT2-FM before DESPOT1
  error('Cannot compute DESPOT2-FM before DESPOT1-HIFI. Do run_despot1() first.');
end

% Try to load the mask, if it exists
if isfield(status, 'mask') && status.mask == 1
  % Load the mask
  mask = load_nifti([dir.MASK status.maskname]);
  disp(['Using user-supplied mask.']);
else
  % Threshold above 0
  mask = ssfp_180(:,:,:,1) > 0;
  disp('Using threshold mask.');
end

% Get size of data
dataSize = size(ssfp_180);

% Mask data
ssfp_0 = ssfp_0 .* repmat(mask, [1 1 1 dataSize(4)]);

% Normalize data by maximum voxel value in masked image
maxVal1  = max(ssfp_0(:));
maxVal2  = max(ssfp_180(:));
maxVal   = max([maxVal1 maxVal2]);

ssfp_0   = ssfp_0   ./ maxVal;
ssfp_180 = ssfp_180 ./ maxVal;

status.despot2_signalScale = maxVal;
disp(['Normalizing SSFP Data By Scale Factor: ' num2str(status.despot2_signalScale)]);

%% EXT B0
if ideal_flag == 1
  disp(['Using ' dir.DESPOT1 'IDEAL-Omega.nii for external B0 map.']);
  opts.ext_omega = img.omega;
  img.omega = [];
end

%% Run DESPOT2-FM for R2/B0 Map
warning off %#ok<WNOFF>
diary('off');
opts.debug = 0;
[pd r2 omega rnrm] = despot2fm_model_fit(reshape(ssfp_0, [dataSize(1)*dataSize(2)*dataSize(3) dataSize(4)]), reshape(ssfp_180, [dataSize(1)*dataSize(2)*dataSize(3) dataSize(4)]), alpha_ssfp, tr_ssfp, r1(:), pd(:), fam(:), opts);
warning on %#ok<WNON>
diary('_mcdespot_log.txt');

% Rehape Results
r2  =   reshape(r2,    [dataSize(1) dataSize(2) dataSize(3)]);
pd  =   reshape(pd,    [dataSize(1) dataSize(2) dataSize(3)]);
omega = reshape(omega, [dataSize(1) dataSize(2) dataSize(3)]);
rnrm  = reshape(rnrm,  [dataSize(1) dataSize(2) dataSize(3)]);

% Save NIfTI
img_nifti_to_nifti(iminv(r2), [dir.BASE dir.SPGR 'spgr_01.nii'], [dir.DESPOT1 'DESPOT2-T2']);
img_nifti_to_nifti(omega,     [dir.BASE dir.SPGR 'spgr_01.nii'], [dir.DESPOT1 'DESPOT2-Omega']);
img_nifti_to_nifti(pd,        [dir.BASE dir.SPGR 'spgr_01.nii'], [dir.DESPOT1 'DESPOT2-PD']);
img_nifti_to_nifti(rnrm,      [dir.BASE dir.SPGR 'spgr_01.nii'], [dir.DESPOT1 'DESPOT2-Rnrm']);

% Plot the center slice
centerSlice = round(dataSize(3) / 2);

% T2 Plot
subplot(1,2,1);
imagesc(iminv(r2(:,:, centerSlice)), [0 1]);
colorbar;
colormap gray;
axis image;
axis off;
title 'T2 DESPOT2 [s]';

% FAM Plot
subplot(1,2,2);
imagesc(omega(:,:, centerSlice), [0 1/tr_ssfp]);
colorbar;
colormap gray;
axis image;
axis off;
title 'Omega (B0 Map) [Hz]';

% Save Figure
% savefig([dir.DESPOT1 'DESPOT2-Result.tif']);

% Done!
time.despot2_end = datetime();
disp(['DESPOT2-FM Complete: ' time.despot2_end]);

% DESPOT2 Flag
status.despot2 = 1;
  
save('_mcdespot_settings', 'tr_spgr', 'tr_irspgr', 'tr_ssfp', 'flags', ...
     'alpha_spgr', 'alpha_afi', 'alpha_irspgr', 'alpha_ssfp', 'ti_irspgr', 'npe_irspgr', ...
     'info_spgr', 'info_irspgr', 'info_afi', 'info_ssfp_0', 'info_ssfp_180', 'info_ideal', 'dir', 'time', 'status');
    
diary('off');

