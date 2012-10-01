% FUNCTION [] = load_mcdespot_vnmr(dirnames)
%
% Load in mcDESPOT data from VARIAN fid data. into
% NIfTI .img/.hdr format.  Naming convention follows
% Sean Deoni's MAC python processing scripts, or Sam's code for
% Matlab processing.
%
% Inputs
%    dirnames - Path to FID files: SPGR, SSFP-180, SSFP-0
%
% Samuel A. Hurley
% University of Wisconsin
% v3.1 1-Dec-2011
%
% Revision History
%    v3.0  Initial Version based on load_mcdespot_data (Oct-2010)
%    v3.1  Change img_3d_to_nifti to img_vnmr_to_nifti (Dec-2011)

function [] = load_mcdespot_vnmr(dirnames)

% Initiate Diary
diary('_mcdespot_log.txt');

VER = 3.1;
VERDATE = '1-Dec-2011';

% Directory Locations (no IR-SPGR on Vnmr1)
dir.BASE     = './originalData/';
dir.SPGR     = 'spgr/';
dir.IRSPGR   = 'irspgr/';
dir.SSFP_0   = 'ssfp_0/';
dir.SSFP_180 = 'ssfp_180/';

% Display banner
disp('=== cpMCDESPOT - Multicomponent Relaxomtery Analysis ===');
disp('     Data Loading Script (mcdespot_load)'                );
disp('     Samuel A. Hurley      shurley@wisc.edu'             );
disp(['    Version ' num2str(VER, '%01.1f') '            ' VERDATE]);
disp('     FOR USE ONLY AT UNIVERSITY OF WISCONSIN.'           );
disp('========================================================');

% Starting time for loading data
time.load_start = datetime();
disp(['Loading Started: ' time.load_start]);

%% O. Acquisition Settings
% Flip angles are no longer hard-coded, read procpar header.

% No IR-SPGR Dataset

%% I. Directories for Images

% Check if directories were specified
if (exist('dirnames', 'var'))
  sts = 1;
else
  % Read in input directories
  [dirnames sts] = spm_select(3, 'dir', 'Select SPGR, SSFP-180, SSFP-0');
end

% If user aborts the dialog...
if sts ~= 1
  return;
end

% Create new directories
mkdir(dir.BASE);
mkdir([dir.BASE dir.SPGR]);
mkdir([dir.BASE dir.SSFP_0]);
mkdir([dir.BASE dir.SSFP_180]);

%% II. Load in PROCPAR Headers

info_spgr     = load_procpar([dirnames(1,:) 'procpar']);
info_ssfp_180 = load_procpar([dirnames(2,:) 'procpar']);
info_ssfp_0   = load_procpar([dirnames(3,:) 'procpar']);

%% III. Grab Info From PROCPAR Headers

% Get TR
tr_spgr      = info_spgr.tr;                    % In s
tr_ssfp      = info_ssfp_0.tr;                  % In s

% IR-SPGR Variables which are not used for Vnmr
tr_irspgr     = [];        %#ok<NASGU>
ti_irspgr     = [];        %#ok<NASGU>
npe_irspgr    = [];        %#ok<NASGU>
alpha_irspgr  = [];
irspgr_fnames = [];
info_irspgr   = [];

% Get flip angles
alpha_spgr    = info_spgr.flip1;
alpha_ssfp    = info_ssfp_0.flip1;

% Display TR & TI
disp( '========================================================');
disp(['SPGR TR:            ' num2str(tr_spgr*1000)    ' ms']);
disp(['SSFP TR:            ' num2str(tr_ssfp*1000)    ' ms']);

% Check Reciever Gain & Presig Attenuator
disp( '========================================================');
disp(['SPGR     Presig/Gain:  ' info_spgr.presig     '/' num2str(info_spgr.gain,     '%02i')]);
disp(['SSFP-180 Presig/Gain:  ' info_ssfp_180.presig '/' num2str(info_ssfp_180.gain, '%02i')]);
disp(['SSFP-0   Presig/Gain:  ' info_ssfp_0.presig   '/' num2str(info_ssfp_0.gain,   '%02i')]);

% Display Flip Angles
disp(['========================================================']);
disp(['SPGR Flip Angles:   ' num2str(alpha_spgr,    '%02.0f ')    ' degrees']);
disp(['SSFP Flip Angles:   ' num2str(alpha_ssfp,    '%02.0f ')    ' degrees']);

% Display the SSFP phase cycling is correct
disp(['========= Check for correct SSFP Phase Cycling =========']); %#ok<NBRAK>
disp([' SSFP 180: ' info_ssfp_180.comment]);
disp([' SSFP 0:   ' info_ssfp_0.comment]);
disp(['========================================================']); %#ok<NBRAK>

% Check that two SSFP have same TR
if info_ssfp_180.tr ~= tr_ssfp
  error('TR of Two SSFP Series Does Not Match!');
end

if (info_ssfp_0.presig ~= info_ssfp_180.presig)
  error('SSFP''s Have Different Presig Setting (High/Low)');
end

if (info_ssfp_0.gain ~= info_ssfp_180.gain)
  error('SSFP''s Have Different Transmitter Gains');
end

% Check that two SSFP have same flip angle
if sum(info_ssfp_180.flip1 == info_ssfp_0.flip1) ~= numel(info_ssfp_180.flip1)
  error('Two SSFP datasets must be acquired with the same flip angles');
end

% Wait a few seconds for the user to check the settings
countdown(5);


%% III. Recon FID Images, Then Convert to NIFTI

disp('FID->NIfTI reconstruction...');

% SPGR
img_spgr = recon_3d(dirnames(1,:));
img_vnmr_to_nifti(img_spgr, info_spgr, [dir.BASE dir.SPGR 'spgr']);

% SSFP-0
img = recon_3d((dirnames(3,:)));
% Check for matching IMG size
if size(img, 1) ~= size(img_spgr, 1) || size(img, 2) ~= size(img_spgr, 2)
  disp('Re-Size SSFP-0 image to match SPGR.');
  img = img_resize(img, [size(img_spgr, 1) size(img_spgr, 2)]);
end
% Save as NIfTI
img_vnmr_to_nifti(img, info_ssfp_0, [dir.BASE dir.SSFP_0 'ssfp_0']);

% SPGR-180
img = recon_3d(strtrim(dirnames(2,:)));
% Check for matching IMG size
if size(img, 1) ~= size(img_spgr, 1) || size(img, 2) ~= size(img_spgr, 2)
  disp('Re-Size SSFP-180 image to match SPGR.');
  img = img_resize(img, [size(img_spgr, 1) size(img_spgr, 2)]);
end

% Save as NIfTI
img_vnmr_to_nifti(img, info_ssfp_180, [dir.BASE dir.SSFP_180 'ssfp_180']);

%% IV Ending time for loading data
time.load_end = datetime();
disp(['Loading Complete: ' time.load_end]);

%% V Save settings for later mcDESPOT processing

% Struct to hold status of processing
status = struct(); %#ok<NASGU>

save('_mcdespot_settings', 'tr_spgr', 'tr_irspgr', 'tr_ssfp', ...
      'alpha_spgr', 'alpha_irspgr', 'alpha_ssfp', 'ti_irspgr', 'npe_irspgr', ...
      'info_spgr', 'info_irspgr', 'info_ssfp_0', 'info_ssfp_180', 'dir', 'time', 'status');

    
diary('off');
