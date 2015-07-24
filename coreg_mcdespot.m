% FUNCTION [] = coreg_mcdespot(target);
%
% FUNCTION to automate loading & coregistration of mcDESPOT
% data using FLIRT
%
% Inputs:
%    Explicit:
%              target - (optional) an external target image for
%              coregistration
%
%    Implicit: Start in same directory as _mcdespot_settings.mat
%
% Outputs:
%    Explicit: None
%    Implicit: New directory called coregisteredData
%
% Samuel A. Hurley
% v5.1 - 26-Jun-2015
%
% Changelog:
%     v3.0 - Initial Version (using v3.0 to match other mcDESPOT commands)  (Jun-2010)
%     v3.1 - Changed directory structure slightly, keep all formats in
%            NIFTI (June-2010)
%     v3.2 - Added input argument to specify an external target image for
%            coregistration (Oct-2010)
%     v3.3 - Tweaked settings, added omat so that transforms are saved (Jan-2011)
%     v3.4 - Updated to be compatible with afi_flag & ideal_flag options (Feb-2012)
%     v3.5 - Updated to daisy-chain coreg of SPGR as well as SSFP-0 (since 6th or 7th
%            SPGR FA ususally show very little contrast, as opposted to PD or T1-w).
%     v5.0 - Skip 4.0 to make syc release versioning
%            Tweaked coreg ranges to reduce chances of SSFP-0 failing. Added support
%            for Bloch-Siegert B1 flag. (Mar-2014)
%     v5.1 - Remove daisy-chained registrations. Use corratio (FSL default)
%            instead of mutual information for similarity metric.

function coreg_mcdespot(target)

%% Initialization

% Initiate Diary
diary('_mcdespot_log.txt');
clc;

VER = 5.1;
VERDATE = '26-Jun-2015';

% Display banner
disp(['=== cpMCDESPOT - Multicomponent Relaxomtery Analysis ===']); %#ok<*NBRAK>
disp(['     Coreg Script - FLIRT (coreg_mcdespot)'              ]);
disp(['     Samuel A. Hurley      shurley@wisc.edu'             ]);
disp(['     Version ' num2str(VER, '%01.1f') '         ' VERDATE]);
disp(['========================================================']);

% Load in mcdespot settings file
load _mcdespot_settings;

% Starting time for coreg data
time.coreg_start = datetime_stamp();
disp(['Coregistration Started: ' time.coreg_start]);

% Make coreg image directories
dir.COREG = './coregisteredData/';

% Create new directories
mkdir(dir.COREG);
mkdir([dir.COREG dir.SPGR]);
mkdir([dir.COREG dir.IRSPGR]);
mkdir([dir.COREG dir.SSFP_0]);
mkdir([dir.COREG dir.SSFP_180]);
mkdir([dir.COREG dir.EXT_CAL]);

% Check Flags for AFI & IDEAL Scans
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

% Check for External Coregistration Target
if exist('target', 'var')
  disp('Using External Coregistration Target');
  ref = target;
  
else
  % Coregistration Reference for SPGR
  disp('Using 1st SPGR Scan as Coregistration Target');
  ref = [dir.BASE dir.SPGR 'spgr_01'];
end

% FLIRT Additional Options
%   bins 256 is default
%   cost corratio is default
%   search_rx/ry/rz - assumes no more than +/- 2 degree rotation
%   set to -nosearch in SSFP-0 section
opts = [' -datatype float  -searchcost corratio -cost corratio -bins 256 -dof 6 -searchrx -2 2' ...
        ' -searchry -2 2 -searchrz -2 2  -coarsesearch 2 -finesearch 1 -interp sinc'];
      
% If mask already exists, specify as additional option
if isfield(status, 'mask') && status.mask == 1 %#ok<NODEF>
  disp(['Using ' status.maskname ' For Brain Ref Weighting']);
  opts = [opts, ' -refweight ' dir.MASK status.maskname];
end

%% Coreg First SPGR
ii = 1;
  disp(['--Flip Angle ' num2str(alpha_spgr(ii), '%02.0f') '--']);
  % Source File
  in  = [dir.BASE dir.SPGR 'spgr_'  num2str(ii, '%02.0f')];
  % Result Filename
  out = [dir.COREG dir.SPGR 'spgr_' num2str(ii, '%02.0f')];
  
  eval(['!flirt -in ' in ' -ref ' ref ' -out ' out ' -omat ' out '.txt' opts]);
  eval(['!fslchfiletype NIFTI ' out]);


%% Coreg Rest of SPGR
disp('== Coregistration of SPGR ==');
for ii = 2:length(alpha_spgr)
  disp(['--Flip Angle ' num2str(alpha_spgr(ii), '%02.0f') '--']);
  
  % Source File
  in  = [dir.BASE dir.SPGR 'spgr_'  num2str(ii, '%02.0f')];
  % Result Filename
  out = [dir.COREG dir.SPGR 'spgr_' num2str(ii, '%02.0f')];
  
  eval(['!flirt -in ' in ' -ref ' ref ' -out ' out ' -omat ' out '.txt' opts]);
  eval(['!fslchfiletype NIFTI ' out]);
end


%% Coreg The FA Calbiration Map

% Set T1w SPGR as ref for IR-SPGR or AFI, since they have more similar constrast
ref = [dir.COREG dir.SPGR 'spgr_' num2str(length(alpha_spgr), '%02.0f')];

if afi_flag == 0
  % IR-SPGR
  disp('== Coregistration of IR-SPGR ==');
  % Source File
  in  = [dir.BASE dir.IRSPGR 'irspgr'];
  % Result Filename
  out = [dir.COREG dir.IRSPGR 'irspgr'];
  eval(['!flirt -in ' in ' -ref ' ref ' -out ' out ' -omat ' out '.txt' opts]);
  eval(['!fslchfiletype NIFTI ' out]);
elseif afi_flag == 1
  % AFI
  disp('== Coregistration of AFI ==');
  disp('--AFI TR 1--');
  % Source File
  in  = [dir.BASE dir.EXT_CAL 'afi_01'];
  % Result Filename
  out = [dir.COREG dir.EXT_CAL 'afi_01'];
  eval(['!flirt -in ' in ' -ref ' ref ' -out ' out ' -omat ' out '.txt' opts]);
  eval(['!fslchfiletype NIFTI ' out]);
  
  % Apply same xform to AFI_02
  disp('--AFI TR 2--');
  in  = [dir.BASE  dir.EXT_CAL 'afi_02'];
  out = [dir.COREG dir.EXT_CAL 'afi_02'];
  eval(['!flirt -in ' in ' -ref ' ref ' -out ' out ' -applyxfm -init ' dir.COREG dir.EXT_CAL 'afi_01.txt' opts]);
  eval(['!fslchfiletype NIFTI ' out]);
  
elseif afi_flag == 2 % Bloch-Siegert - apply identity XFORM to make resolution match
  disp('== Coreg BS B1 Map ==');
  disp('--BS (Mag)--');
  
  % Apply XFORM to resample BS into same space as SPGR
  ref = [dir.COREG dir.SPGR 'spgr_01'];
  % Source File
  in  = [dir.BASE dir.EXT_CAL 'afi_02'];
  % Result Filename
  out = [dir.COREG dir.EXT_CAL 'afi_02'];
  
  % Write out Identity Xform
  write_xform(eye(4), [dir.COREG dir.EXT_CAL 'XFORM.txt']);
  
  eval(['!flirt -in ' in ' -ref ' ref ' -out ' out ' -applyxfm -init ' dir.COREG dir.EXT_CAL 'XFORM.txt' opts]);
  eval(['!fslchfiletype NIFTI ' out]);
  
  % Apply same xform to AFI_02
  disp('--BS (Map)--');
  in  = [dir.BASE dir.EXT_CAL  'afi_01'];
  out = [dir.COREG dir.EXT_CAL 'afi_01'];
  
  eval(['!flirt -in ' in ' -ref ' ref ' -out ' out ' -applyxfm -init ' dir.COREG dir.EXT_CAL 'XFORM.txt' opts]);
  eval(['!fslchfiletype NIFTI ' out]);

end

%% Coreg IDEAL - Apply Identity XFORM to IDEAL to Make Sure Resolution Matches
if ideal_flag == 1
  disp('== Apply I4 Transform to IDEAL B0 Map ==');
  % Write out Identity Xform
  write_xform(eye(4), [dir.COREG dir.EXT_CAL 'IDEAL.txt']);
  in  = [dir.BASE dir.EXT_CAL 'IDEAL-Omega'];
  out = [dir.COREG dir.EXT_CAL 'IDEAL-Omega'];
  eval(['!flirt -in ' in ' -ref ' ref ' -out ' out ' -applyxfm -init ' dir.COREG dir.EXT_CAL 'IDEAL.txt' opts]);
  eval(['!fslchfiletype NIFTI ' out]);
end

%% Coreg 1st SSFP-180

% Next, coregister the first SSFP-180 image to the PD-weighted SPGR
% Because they have similar contrast, then coreg all SSFP images to the
% first SSFP-180 image.

% Reference Image - PDw SPGR
ref = [dir.COREG dir.SPGR 'spgr_01'];

disp('== Coregistration of SSFP-180 ==');
disp(['--Flip Angle ' num2str(alpha_ssfp(1), '%02.0f') ' Cycle: 180--']);

% Source File
in  = [dir.BASE dir.SSFP_180 'ssfp_180_01'];
% Result Filename
out = [dir.COREG dir.SSFP_180 'ssfp_180_01'];

% Apply init matrix
eval(['!flirt -in ' in ' -ref ' ref ' -out ' out ' -omat ' out '.txt' opts]);
eval(['!fslchfiletype NIFTI ' out]);

%% Coreg Rest of SSFP-180

% Update Reference Image
ref = [dir.COREG dir.SSFP_180 'ssfp_180_01'];

for ii = 2:length(alpha_ssfp)
  disp(['--Flip Angle ' num2str(alpha_ssfp(ii), '%02.0f') ' Cycle: 180--']);
  
  % Source File
  in  = [dir.BASE  dir.SSFP_180 'ssfp_180_' num2str(ii, '%02.0f')];
  % Result Filename
  out = [dir.COREG dir.SSFP_180 'ssfp_180_' num2str(ii, '%02.0f')];
  
  eval(['!flirt -in ' in ' -ref ' ref ' -out ' out ' -omat ' out '.txt' opts]);
  eval(['!fslchfiletype NIFTI  ' out]);
end

%% Coreg 1st SSFP-0

disp('== Coregistration of SSFP-0 Data ==');
ii = 1;

% Update Reference Image & Inital Transform to 1st SSFP-180 image
% Add -nosearch option for SSFP-0 Coreg to prevent large rotations/shifts in banded images
ref   = [dir.COREG dir.SSFP_180 'ssfp_180_01'];
xfm   = [dir.COREG dir.SSFP_180 'ssfp_180_01'];
opts1 = [opts ' -nosearch -init ' xfm '.txt '];

disp(['--Flip Angle ' num2str(alpha_ssfp(1), '%02.0f') ' Cycle: 0--']);

% Source File
in  = [dir.BASE  dir.SSFP_0 'ssfp_0_' num2str(ii, '%02.0f')];
% Result Filename
out = [dir.COREG dir.SSFP_0 'ssfp_0_' num2str(ii, '%02.0f')];

% Apply init matrix
eval(['!flirt -in ' in ' -ref ' ref ' -out ' out ' -omat ' out '.txt' opts1]);
eval(['!fslchfiletype NIFTI  ' out  ]);

%% Coreg Rest of SSFP-0

% Update Reference Image & Initial Transform to 1st SSFP-0
% Add -nosearch option for SSFP-0 Coreg to prevent large rotations/shifts in low SNR images
ref   = [dir.COREG dir.SSFP_0 'ssfp_0_' num2str((ii), '%02.0f')];
xfm   = [dir.COREG dir.SSFP_0 'ssfp_0_' num2str((ii), '%02.0f')];
opts1 = [opts ' -nosearch -init ' xfm '.txt '];

for ii = 2:length(alpha_ssfp)
  
  disp(['--Flip Angle ' num2str(alpha_ssfp(ii), '%02.0f') ' Cycle: 0--']);
  
  % Source File
  in  = [dir.BASE  dir.SSFP_0 'ssfp_0_' num2str(ii, '%02.0f')];
  % Result Filename
  out = [dir.COREG dir.SSFP_0 'ssfp_0_' num2str(ii, '%02.0f')];

  eval(['!flirt -in ' in ' -ref ' ref ' -out ' out ' -omat ' out '.txt' opts1]);
  eval(['!fslchfiletype NIFTI  ' out]);
end

%% Finish Up

time.coreg_end = datetime_stamp();
disp(['Coreg Complete: ' time.coreg_end]);

% Coreg Flag
status.coreg = 1;

save('_mcdespot_settings', 'tr_spgr', 'tr_irspgr', 'tr_ssfp', 'flags', ...
     'alpha_spgr', 'alpha_afi', 'alpha_irspgr', 'alpha_ssfp', 'ti_irspgr', 'npe_irspgr', ...
     'info_spgr', 'info_irspgr', 'info_afi', 'info_ssfp_0', 'info_ssfp_180', 'info_ideal', 'dir', 'time', 'status');
    
diary('off');