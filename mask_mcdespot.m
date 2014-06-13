% FUNCTION [] = mask_mcdespot(flip, prog);
%
% FUNCTION to automate brain extraction
%
% Inputs:
%    Explicit: prog - 0=fsl's BET, 1=AFNI's 3dSkullStrip
%              flip - which SPGR flip angle to use
%    Implicit: Start in same directory as _mcdespot_settings.mat
%
% Outputs:
%    Explicit: None
%    Implicit: New directory called maskData
%
% Samuel A. Hurley
% v3.3 28-Feb-2012
%
% Changelog:
%     v3.0 - Initial Version (using v3.0 to match other mcDESPOT commands)
%     v3.1 - Options to choose which SPGR image (flip angle) to use  (May-2010)
%     v3.2 - Added additional arguemnts to BET (Oct-2010)
%     v3.3 - Updated to be compatiable with afi_flag & ideal_flag options  (Feb-2012)
%     v5.0 - Updated to use gzipped mask (instead of .nii mask) - easier for FSLView (Jun-2014)

function mask_mcdespot(flip, prog)

% Constants -- BET Tweaks
THRESH_OM =  0.5;
THRESH_G  = -0.1;

% Initiate Diary
diary('_mcdespot_log.txt');

% Default options
if ~exist('flip', 'var')
  % Default is PD-Weighted SPGR
  flip = 1;
end
if ~exist('prog', 'var')
  % Default is FSL's BET tool
  prog = 0;
end

% Display Banner
VER     = 5.0;
VERDATE = 'Jun-2014';

% Display banner
disp(['=== cpMCDESPOT - Multicomponent Relaxomtery Analysis ===']); %#ok<*NBRAK>
disp(['     Brain Extraction Script (mask_mcdespot)'            ]);
disp(['     Samuel A. Hurley      shurley@wisc.edu'             ]);
disp(['     Pouria Mossahebi      mossahebi@wisc.edu'           ]);
disp(['     Version ' num2str(VER, '%01.1f') '         ' VERDATE]);
disp(['     FOR USE ONLY AT UNIVERSITY OF WISCONSIN.'           ]);
disp(['========================================================']);

% Load in mcdespot settings file
load _mcdespot_settings;

% Starting time for coreg data
time.mask_start = datetime();
disp(['Masking/BET Started: ' time.mask_start]);

% Define mask location
dir.MASK = './maskData/';
mkdir(dir.MASK);

% Check correct dir
if isfield(status, 'coreg') && status.coreg == 1 %#ok<NODEF>
  basedir = dir.COREG;
else
  basedir = dir.BASE;
end

% Masking Reference for BET or 3dSkullStrip
ref = [basedir dir.SPGR 'spgr_' num2str(flip, '%02.0f')];

flipstr = num2str(alpha_spgr(flip), '%02.0f');

% Brain extract the dataset
disp('== Brain Extracting SPGR ==');
if prog == 0
  disp('Using FSL BET...');
  eval(['!bet2 ' ref ' ' fullfile(dir.MASK, ['bet_fa' flipstr]) ' -m -r 90 -n -f ' num2str(THRESH_OM)  ' -g ' num2str(THRESH_G) ' -t']);
  %V5.0 - Don't unzip FSL mask eval(['!gunzip ' fullfile(dir.MASK, ['bet_fa' flipstr '_mask.nii.gz'])]);
  status.maskname = ['bet_fa' flipstr '_mask.nii.gz'];

elseif prog == 1
  disp('Using AFNI 3dSkullStrip...');
  eval(['!3dSkullStrip -input ' ref '.nii -prefix ' fullfile(dir.MASK, ['afni_fa' flipstr '_mask']) ' -mask_vol']);
  eval(['!3dAFNItoNIFTI ' fullfile(dir.MASK, ['afni_fa' flipstr '_mask+orig'])]);
  eval(['!mv afni_fa' flipstr '_mask.nii ' dir.MASK 'afni_fa' flipstr '_mask.nii']);
  eval(['!rm -rf ' dir.MASK 'afni_fa' flipstr '_mask+orig*']);
  status.maskname = ['afni_fa' flipstr '_mask.nii'];
  
else
  error('Incorrect prog argument specified');
end

% Load the mask back in and save a screenshot
mask = load_nifti([dir.MASK status.maskname]);
spgr = load_nifti([basedir dir.SPGR 'spgr_01.nii']);

dataSize = size(spgr);
centerSlice = round(dataSize ./ 2);

subplot(1,3,1);
cs = centerSlice(1);
imagesc(squeeze(spgr(cs,:,:)).*(squeeze(mask(cs,:,:)) + ~squeeze(mask(cs,:,:))*.2))
colormap gray;
axis image;
axis off;

subplot(1,3,2);
cs = centerSlice(2);
imagesc(squeeze(spgr(:,cs,:)).*(squeeze(mask(:,cs,:)) + ~squeeze(mask(:,cs,:))*.2))
colormap gray;
axis image;
axis off;

subplot(1,3,3);
cs = centerSlice(3);
imagesc(spgr(:,:,cs).*(mask(:,:,cs) + ~mask(:,:,cs)*.2))
colormap gray;
axis image;
axis off;

saveas(gcf, [dir.MASK 'Mask.tif']);

% Done!
time.mask_end = datetime();
disp(['Masking/BET Complete: ' time.mask_end]);

% Coreg Flag
status.mask = 1;

% Save the program used & the SPGR image used
% to generate the mask
status.masktype  = prog;
status.maskimg   = flip;
  
save('_mcdespot_settings', 'tr_spgr', 'tr_irspgr', 'tr_ssfp', 'flags', ...
     'alpha_spgr', 'alpha_afi', 'alpha_irspgr', 'alpha_ssfp', 'ti_irspgr', 'npe_irspgr', ...
     'info_spgr', 'info_irspgr', 'info_afi', 'info_ssfp_0', 'info_ssfp_180', 'info_ideal', 'dir', 'time', 'status');
  
diary('off');
