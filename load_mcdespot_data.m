% FUNCTION [] = load_mcdespot_data(dirnames, flipset, flags)
%
% Load in mcDESPOT data from DICOM files, convert to
% ANALYZE .img/.hdr format.  Naming convention follows
% Sean Deoni's MAC python processing scripts, or Sam's code for
% Matlab processing.
%
% Inputs
%    dirnames - Path to DICOM files: SPGR, IR-SPGR/AFI, SSFP-180, SSFP-0
%    flipset  - Specify flip angle set 1 = Deoni 2 = Ollinger
%    flags
%             - .afi: use AFI instaed of DESPOT1-HIFI for FA mapping (0=HIFI 1=AFI 2=BS)
%             - .ideal: use IDEAL instead of DESPOT2-FM for B0 mapping
%             - .vnmr: load VnmrJ pre-clinical data instead of GE DICOM
%             - .reorient: reorient all image series to Axial (improves  BET)
%
% Samuel A. Hurley
% University of Wisconsin
% v3.9 12-Sept-2012
%
% Revision History
%    Based on load_mcdespot_data (for waisman data)
%    v1.0 to read data saved from WIMR MR scanners (I***.dcm format) (2-Dec-2010)
%    v2.0 Waisman & WIMR Loading is now incorporated into same file
%         Now saves settings of IR-SPGR as well.
%         Re-organized file into cell arrays
%    v2.1 Options to specify directories @ command line
%         and auto-detect file format (WIMR/WAISMAN)
%    v2.2 Converts ANALYZE files to 32-bit float, to be compatible with
%         Sean's MAC codes (Feb-2010)
%    v2.3 Changed 32bit FLOAT conversion to use SPM
%    v3.0 Removed Deoni's File Naming Convention for Something simplier
%         Now saves entire DICOM header in the _mcdespot_settings folder
%         Removed all calls to progressbar() (because it steals window focus)
%         (May-2010)
%    v3.1 Replaced spm_dicom_convert with img_dcm_to_nifti conversion tool
%         Flip-flop SSFP-0 and SSFP-180 in load dialog (Jun-2010).
%    v3.2 Transmitter gain check (Jun-2010)
%    v3.3 Added support for Waisman all-in-one zip files (WAISMAN2) (Jul-2010)
%    v3.4 Removed ftype (since it is detected automatically)  Added flag
%         for different flip angle sets (so we can use Jonh's optimized set
%         as well) (Jan-2011)
%    v3.5 Removed resize data for SSFP-180 ZIP2 stuff (because this breaks the new
%         img_dcm_to_nifti command version (April-2011)
%    v3.6 Creates 4D NIfTI Files for the sake of brevity (Jul-2011)
%    v3.7 Reorients to Axial (makes BET work, and easier to view in FSLView) (13-Feb-2012)
%    v3.8 Combined load_mcdespot_vnmr with load_mcdespot_data (to cover AFI vs hifi case) (27-Feb-2012)
%    v3.9 Uses check_gain to handle if SSFP phase cycles end up with different reciever gain settings (12-Sept-2012)
%
%    v5.0 NOTE: Jump to 5.0 to keep source maj rev in sync.
%               Added support for Bloch-Siegert (BS) B1 correction (in addition to AFI and HIFI)

function [] = load_mcdespot_data(dirnames, flipset, flags)

% Initiate Diary
diary('_mcdespot_log.txt');
%clc;

VER = 5.0;
VERDATE = '3-Mar-2014';

% Directory Locations
dir.BASE     = './originalData/';
dir.SPGR     = 'spgr/';
dir.IRSPGR   = 'irspgr/';
dir.SSFP_0   = 'ssfp_0/';
dir.SSFP_180 = 'ssfp_180/';
dir.EXT_CAL  = 'calibrationMaps/';

% Display banner
disp(['=== cpMCDESPOT - Multicomponent Relaxomtery Analysis ===']); %#ok<*NBRAK>
disp(['     Data Loading Script (mcdespot_load)'                ]);
disp(['     Samuel A. Hurley      shurley@wisc.edu'             ]);
disp(['     Pouria Mossahebi      mossahebi@wisc.edu'           ]);
disp(['     Version ' num2str(VER, '%01.1f') '         ' VERDATE]);
disp(['     FOR USE ONLY AT UNIVERSITY OF WISCONSIN.'           ]);
disp(['========================================================']);

fprintf('Processed by: ' );
!whoami
fprintf('On host:      ');
!hostname

% Starting time for loading data
time.load_start = datetime();
disp(['Loading Started: ' time.load_start]);

%% O. Acquisition Settings
% # Flip angle scale factors are hard-coded into DESPOT PSD. Search for (multiflip_flag == 3)
if ~check_var('flipset')
  flipset = 1; % Default is Deoni
  disp('FA Set: Defaulting to Deoni 8FA');
end

if flipset == 1
  flips_spgr = [0.16666667 0.22222222 0.27777778 0.33333333 0.38888889 0.50000000 0.72222222 1.00];
  flips_ssfp = [0.17647059 0.23529412 0.30882353 0.39705882 0.48529411 0.58823529 0.75000000 1.00];
  disp('FA SET: Deoni 8FA');
elseif flipset == 2
  flips_spgr = [0.16666667 0.22222222 0.27777778 0.33333333 0.38888889 0.50000000 0.72222222 1.00];
  flips_ssfp = [0.04       0.10       0.20       0.30       0.40       0.60       0.80       1.00];
  disp('FANG Knee FAs');
else
  error('Additional flip angle sets are not supported');
end

% # of phase encodes per IR for DESPOT-HiFi -- Can Read from DICOM header
% npe_irspgr = 96; %#ok<NASGU>
% disp('IR-SPGR: Default to 96PE per inversion');

% Check Flags Struct
if ~check_var('flags')
  flags = struct();
end

% Check Flags for AFI, IDEAL, VNMR, and Reorient
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

if ~isfield(flags, 'vnmr')
  vnmr_flag = 0;
else
  vnmr_flag = flags.vnmr;
end

if ~isfield(flags, 'reorient')
  reorient_flag = 0;
else
  reorient_flag = flags.reorient;
end



%% I. Directories for Images

% Check if directories were specified
if check_var('dirnames')
  sts = 1;
else
  % Read in input directories
  if afi_flag == 1 && ideal_flag == 1
    [dirnames sts] = spm_select(5, 'dir', 'Select SPGR, AFI, SSFP-180, SSFP-0, IDEAL');
  elseif afi_flag == 2 && ideal_flag == 1
    [dirnames sts] = spm_select(5, 'dir', 'Select SPGR, Bloch-Siegert, SSFP-180, SSFP-0, IDEAL');
  elseif afi_flag == 1 && ideal_flag == 0
    [dirnames sts] = spm_select(4, 'dir', 'Select SPGR, AFI, SSFP-180, SSFP-0');
  elseif afi_flag == 2 && ideal_flag == 0
    [dirnames sts] = spm_select(4, 'dir', 'Select SPGR, Bloch-Siegert, SSFP-180, SSFP-0');
  elseif ideal_flag == 1
    [dirnames sts] = spm_select(5, 'dir', 'Select SPGR, IR-SPGR, SSFP-180, SSFP-0, IDEAL');
  else
    [dirnames sts] = spm_select(4, 'dir', 'Select SPGR, IR-SPGR, SSFP-180, SSFP-0');
  end
end

if sts ~= 1
  disp(['Loading aborted at ' datetime]);
  diary('off');
  return;
end

% Convert dirnames to cell array
dirnames1 = dirnames;
dirnames  = cell(0);

for ii = 1:size(dirnames1, 1)
  dirnames{ii} = strtrim(dirnames1(ii,:));
end

% Create new directories
mkdir(dir.BASE);
mkdir([dir.BASE dir.SPGR]);
mkdir([dir.BASE dir.SSFP_0]);
mkdir([dir.BASE dir.SSFP_180]);

% Create External Cal Directory
if afi_flag > 0 || ideal_flag == 1
  mkdir([dir.BASE dir.EXT_CAL]);
end

% Create IRSPGR Directory
if afi_flag == 0
  mkdir([dir.BASE dir.IRSPGR]);
end

%% II. Load in DICOM Images & Headers
[img_spgr     info_spgr]     = load_dicom(dirnames{1});
[img_ssfp_180 info_ssfp_180] = load_dicom(dirnames{3});
[img_ssfp_0   info_ssfp_0]   = load_dicom(dirnames{4});

if afi_flag == 0
  [img_irspgr info_irspgr]   = load_dicom(dirnames{2});
  info_afi    = [];
elseif afi_flag == 1
  [img_afi    info_afi]      = load_dicom(dirnames{2},[],1); % AFI Flag
  info_irspgr = [];
elseif afi_flag == 2 % SAH 5.0 case for Bloch-Siegert. Keep 'afi' name
  [img_afi    info_afi]      = load_dicom(dirnames{2},[],0);
  info_irspgr = [];
end

if ideal_flag == 1
  [img_ideal  info_ideal]    = load_dicom(dirnames{5});
else
  info_ideal = [];
end

%% III. Grab Info From DICOM Headers

% Get TR
tr_spgr    = info_spgr.RepetitionTime;        % In ms
tr_ssfp    = info_ssfp_0.RepetitionTime;      % In ms

% Get flip angles
alpha_spgr    = info_spgr.FlipAngle .* flips_spgr;
alpha_ssfp    = info_ssfp_0.FlipAngle .* flips_ssfp;

% IR-SPGR Parameters
if afi_flag == 0
  tr_irspgr  = info_irspgr.RepetitionTime;    % In ms;
  alpha_irspgr  = info_irspgr.FlipAngle;
  
  alpha_afi = [];
  
  % Get TI & NPE for the IR-SPGR
  ti_irspgr  = info_irspgr.InversionTime;  % In MS
  npe_irspgr = double(round(info_irspgr.AcquisitionMatrix(3).*info_irspgr.PercentPhaseFieldOfView/200));
else
  alpha_afi = info_afi.FlipAngle;
  
  tr_irspgr    = [];
  alpha_irspgr = [];
  ti_irspgr    = [];
  npe_irspgr   = [];
end

% Display TR & TI
disp( '========================================================');
disp(['SPGR TR:            ' num2str(tr_spgr)    ' ms']);
disp(['SSFP TR:            ' num2str(tr_ssfp)    ' ms']);
if afi_flag == 0
  disp(['IR-SPGR TR:         ' num2str(tr_irspgr)  ' ms']);
  disp(['IR-SPGR TI:         ' num2str(ti_irspgr)  ' ms']);
elseif afi_flag == 1
  disp(['AFI TR:             ' num2str(info_afi.RepetitionTime) ' ms']);
else
  disp(['BS  TR:             ' num2str(info_afi.RepetitionTime) ' ms']);
end

% Check Reciever/Transmitter Gains
disp( '========================================================');
disp(['SPGR     R1/R2/TG:  ' num2str(info_spgr.Private_0019_108a,     '%02i') '/' num2str(info_spgr.Private_0019_108b,     '%02i') '/' num2str(info_spgr.Private_0019_1094,     '%02i')]);
if afi_flag == 0
  disp(['IRSPGR   R1/R2/TG:  ' num2str(info_irspgr.Private_0019_108a, '%02i') '/' num2str(info_irspgr.Private_0019_108b,   '%02i') '/' num2str(info_irspgr.Private_0019_1094,   '%02i')]);
elseif afi_flag == 1
  disp(['AFI      R1/R2/TG:  ' num2str(info_afi.Private_0019_108a,    '%02i') '/' num2str(info_afi.Private_0019_108b,      '%02i') '/' num2str(info_afi.Private_0019_1094,      '%02i')]);
else
  disp(['BS       R1/R2/TG:  ' num2str(info_afi.Private_0019_108a,    '%02i') '/' num2str(info_afi.Private_0019_108b,      '%02i') '/' num2str(info_afi.Private_0019_1094,      '%02i')]);
end
disp(['SSFP-180 R1/R2/TG:  ' num2str(info_ssfp_180.Private_0019_108a, '%02i') '/' num2str(info_ssfp_180.Private_0019_108b, '%02i') '/' num2str(info_ssfp_180.Private_0019_1094, '%02i')]);
disp(['SSFP-0   R1/R2/TG:  ' num2str(info_ssfp_0.Private_0019_108a,   '%02i') '/' num2str(info_ssfp_0.Private_0019_108b,   '%02i') '/' num2str(info_ssfp_0.Private_0019_1094,   '%02i')]);

% Display Flip Angles
disp(['========================================================']);
disp(['SPGR Flip Angles:   ' num2str(alpha_spgr,     '%02.0f ')    ' degrees']);
disp(['SSFP Flip Angles:   ' num2str(alpha_ssfp,     '%02.0f ')    ' degrees']);
if afi_flag == 0
  disp(['IR-SPGR Flip Angle: ' num2str(alpha_irspgr, '%02.0f ')    ' degrees']);
  disp(['IR-SPGR NumPhaseEnc:' num2str(npe_irspgr,   '%03.0f ')              ]);
elseif afi_flag == 1
  disp(['AFI Flip Angle:     ' num2str(alpha_afi,    '%02.0f') ' degrees']);
else
  disp(['BS  Flip Angle:     ' num2str(alpha_afi,    '%02.0f') ' degrees']);
end

% Display the SSFP phase cycling is correct
disp(['========= Check for correct SSFP Phase Cycling =========']); %#ok<NBRAK>
disp([' SSFP 180: ' info_ssfp_180.SeriesDescription]);
disp([' SSFP 0:   ' info_ssfp_0.SeriesDescription]);
disp(['========================================================']); %#ok<NBRAK>

% Check that two SSFP have same TR
if info_ssfp_0.RepetitionTime ~= info_ssfp_180.RepetitionTime
  warning('Two SSFP datasets must be acquired with same TR.');
end

if afi_flag == 0
  % Check that IR was setup for this sequence
  if ti_irspgr == 0
    warning('IR was disabled for IR-SPGR, or the SPGR and IR-SPGR sequences are mixed up.');
  end
  
  % Check that SPGR-IRSPGR and two SSFP have same R1/R2/TG
  if (info_spgr.Private_0019_108a ~= info_irspgr.Private_0019_108a) || (info_spgr.Private_0019_108b ~= info_irspgr.Private_0019_108b)
  error('SPGR and IR-SPGR Have Different Reciever Gains');
  end

  if (info_spgr.Private_0019_1094 ~= info_irspgr.Private_0019_1094)
    disp('SPGR and IR-SPGR Have Different Transmitter Gains');
  end
elseif afi_flag == 1
  % Check that AFI FA is 55
  if info_afi.FlipAngle ~= 55
    warning('AFI Flip Angle Should Be 55 degrees');
  end
  
  % Check that SPGR & AFI Have Same TG
  if (info_spgr.Private_0019_1094 ~= info_afi.Private_0019_1094)
    warning('SPGR and AFI do not have same transmit gain');
  end
end
  
% Check That SSFP Reciever & Transmitter Gains Match
if (info_ssfp_0.Private_0019_108a ~= info_ssfp_180.Private_0019_108a) || (info_ssfp_0.Private_0019_108b ~= info_ssfp_180.Private_0019_108b)
  warning('SSFP''s Have Different Reciever Gains. Attempting to correct');
  % Make ssfp-0 match gain of ssfp-180
  img_ssfp_0 = check_gain(img_ssfp_0, info_ssfp_0, info_ssfp_180);
end

if (info_ssfp_0.Private_0019_1094 ~= info_ssfp_180.Private_0019_1094)
  warning('SSFP''s Have Different Transmitter Gains');
end

% Check that two SSFP have same flip angle
if info_ssfp_0.FlipAngle ~= info_ssfp_180.FlipAngle
  error('Two SSFP datasets must be acquired with same flip angles.');
end

% Wait a few seconds for the user to check the settings
countdown(5);


%% III. Convert images into .NII

disp('DICOM->NIfTI conversion...');

% Call img_dcm_to_nifti to convert to a .nii file
% Note that this command automatically converts to 32-BIT FLOAT,
% so it is compatible with Sean's DESPOT-HIFI precompiled binaries

% Also note that load_dicom automatically sorts the series into a 4-D
% matrix based on flip angle, and automatically un-zips and re-zips Waisman
% Files

% DICOM Version

if vnmr_flag == 0
  
  % SPGR
  img_dcm_to_nifti(img_spgr,     info_spgr,     [dir.BASE dir.SPGR     'spgr'],     0);
  
  % IR-SPGR or AFI
  if afi_flag == 0
    img_dcm_to_nifti(img_irspgr, info_irspgr,   [dir.BASE dir.IRSPGR   'irspgr'],   0);
  else
    img_dcm_to_nifti(img_afi,    info_afi,      [dir.BASE dir.EXT_CAL  'afi'   ],   0);
  end
  
  % SSFP-180
  img_dcm_to_nifti(img_ssfp_180, info_ssfp_180, [dir.BASE dir.SSFP_180 'ssfp_180'], 0);
  
  % SSFP-0
  img_dcm_to_nifti(img_ssfp_0,   info_ssfp_0,   [dir.BASE dir.SSFP_0   'ssfp_0'],   0);
  
  % IDEAL
  if ideal_flag == 1
    img_dcm_to_nifti(img_ideal,  info_ideal,    [dir.BASE dir.EXT_CAL  'IDEAL-Omega'], 0);
  end
  
elseif vnmr_flag == 1
  
  disp('FID->NIfTI reconstruction...');
  
  % SPGR
  img_spgr = recon_3d(dirnames(1,:));
  img_vnmr_to_nifti(img_spgr, info_spgr, [dir.BASE dir.SPGR 'spgr']);
  
  % IR-SPGR or AFI
  if afi_flag == 0
    img_irspgr = recon_3d(dirnames(2,:));
    img_vnmr_to_nifti(img_irspgr, info_irspgr, [dir.BASE dir.IRSPGR 'irspgr']);
  else
    img_afi = recon_3d(dirnames(2,:));
    fam     = afi_standard(img_afi(:,:,:,1,1), img_afi(:,:,:,1,2), 5);
    img_vnmr_to_nifti(fam,    info_afi,    [dir.BASE dir.EXT_CAL, 'AFI-FAM']);
  end
  
  % SPGR-180
  img = recon_3d(strtrim(dirnames(3,:)));
  % Check for matching IMG size
  if size(img, 1) ~= size(img_spgr, 1) || size(img, 2) ~= size(img_spgr, 2)
    disp('Re-Size SSFP-180 image to match SPGR.');
    img = img_resize(img, [size(img_spgr, 1) size(img_spgr, 2)]);
  end
  
  % SSFP-0
  img = recon_3d((dirnames(4,:)));
  % Check for matching IMG size
  if size(img, 1) ~= size(img_spgr, 1) || size(img, 2) ~= size(img_spgr, 2)
    disp('Re-Size SSFP-0 image to match SPGR.');
    img = img_resize(img, [size(img_spgr, 1) size(img_spgr, 2)]);
  end
  
  % Save as NIfTI
  img_vnmr_to_nifti(img, info_ssfp_180, [dir.BASE dir.SSFP_180 'ssfp_180']);
  img_vnmr_to_nifti(img, info_ssfp_0, [dir.BASE dir.SSFP_0 'ssfp_0']);
  
end

%% IV. Reorient All Images to std
%  So that they view correctly on FSLView, and so that BET2 works properly

if reorient_flag == 1
  
  fprintf('Reorient SPGR:');
  for ii = 1:length(alpha_spgr)
    progressbar(ii/length(alpha_spgr));
    invol = [dir.BASE dir.SPGR 'spgr_' num2str(ii, '%02.0f') '.nii'];
    eval(['!fslswapdim ' invol ' RL PA IS ' invol]);
    eval(['!rm -f ' invol]);
    eval(['!gunzip ' invol '.gz']);
  end
  
  fprintf('Reorient SSFP:');
  for ii = 1:length(alpha_ssfp)
    progressbar(ii/length(alpha_ssfp));
    invol = [dir.BASE dir.SSFP_0   'ssfp_0_'   num2str(ii, '%02.0f') '.nii'];
    eval(['!fslswapdim ' invol ' RL PA IS ' invol]);
    eval(['!rm -f ' invol]);
    eval(['!gunzip ' invol '.gz']);
    invol = [dir.BASE dir.SSFP_180 'ssfp_180_' num2str(ii, '%02.0f') '.nii'];
    eval(['!fslswapdim ' invol ' RL PA IS ' invol]);
    eval(['!rm -f ' invol]);
    eval(['!gunzip ' invol '.gz']);
  end
  
  if afi_flag == 0
    fprintf('Reorient IR-SPGR:');
    progressbar(0);
    invol = [dir.BASE dir.IRSPGR 'irspgr.nii'];
    eval(['!fslswapdim ' invol ' RL PA IS ' invol]);
    eval(['!rm -f ' invol]);
    eval(['!gunzip ' invol '.gz']);
    progressbar(1);
  elseif afi_flag == 1
    fprintf('Reorient AFI:');
    progressbar(0);
    invol = [dir.BASE dir.EXT_CAL 'afi_01.nii'];
    eval(['!fslswapdim ' invol ' RL PA IS ' invol]);
    eval(['!rm -f ' invol]);
    eval(['!gunzip ' invol '.gz']);
    
    progressbar(.5);
    invol = [dir.BASE dir.EXT_CAL 'afi_02.nii'];
    eval(['!fslswapdim ' invol ' RL PA IS ' invol]);
    eval(['!rm -f ' invol]);
    eval(['!gunzip ' invol '.gz']);
    progressbar(1);
  end
  
  % SAH 5.0: Note BS should always be acquired in axial, so no need to reorient
  
  if ideal_flag == 1
    fprintf('Reorient IDEAL: ');
    progressbar(0);
    invol = [dir.BASE dir.EXT_CAL 'IDEAL-Omega.nii'];
    eval(['!fslswapdim ' invol ' RL PA IS ' invol]);
    eval(['!rm -f ' invol]);
    eval(['!gunzip ' invol '.gz']);
    progressbar(1);
  end
  
end


%% V Ending time for loading data
time.load_end = datetime();
disp(['Loading Complete: ' time.load_end]);

%% VI Save settings for later mcDESPOT processing
% Convert to seconds for Sam's algorithm
tr_spgr   = tr_spgr   / 1000;     %#ok<NASGU>
tr_irspgr = tr_irspgr / 1000;     %#ok<NASGU>
tr_ssfp   = tr_ssfp   / 1000;     %#ok<NASGU>

% Struct to hold status of processing
status = struct(); %#ok<NASGU>

save('_mcdespot_settings', 'tr_spgr', 'tr_irspgr', 'tr_ssfp', 'flags', ...
     'alpha_spgr', 'alpha_afi', 'alpha_irspgr', 'alpha_ssfp', 'ti_irspgr', 'npe_irspgr', ...
     'info_spgr', 'info_irspgr', 'info_afi', 'info_ssfp_0', 'info_ssfp_180', 'info_ideal', 'dir', 'time', 'status');
    
diary('off');
