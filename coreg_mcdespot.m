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
% v3.4 28-Feb-2012
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

function coreg_mcdespot()

% Initiate Diary
diary('_mcdespot_log.txt');
clc;

VER = 3.4;
VERDATE = '28-Feb-2012';

% Display banner
disp(['=== cpMCDESPOT - Multicomponent Relaxomtery Analysis ===']); %#ok<*NBRAK>
disp(['     Coreg Script - FLIRT (coreg_mcdespot)'              ]);
disp(['     Samuel A. Hurley      shurley@wisc.edu'             ]);
disp(['     Pouria Mossahebi      mossahebi@wisc.edu'           ]);
disp(['     Version ' num2str(VER, '%01.1f') '         ' VERDATE]);
disp(['     FOR USE ONLY AT UNIVERSITY OF WISCONSIN.'           ]);
disp(['========================================================']);

% Load in mcdespot settings file
load _mcdespot_settings;

% Starting time for coreg data
time.coreg_start = datetime();
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
%   searchx/y/z - assumes no more than +/- 5 degree shift
opts = [' -datatype float  -searchcost mutualinfo -cost mutualinfo -bins 256 -dof 6 -searchrx -6 6' ...
        ' -searchry -6 6 -searchrz -6 6  -coarsesearch 2 -finesearch 1 -interp sinc'];
      
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
  
  eval(['!flirt -in ' in ' -ref ' ref ' -out ' out ' -omat ' out '.txt' opts1]);
  eval(['!fslchfiletype NIFTI ' out]);



%% Coreg Rest of SPGR
disp('== Coregistration of SPGR ==');
for ii = 2:length(alpha_spgr)
  
  % Daisy-chain coreg of high FA scans -- Hack-ey way, fix later
  ref   = [dir.COREG dir.SPGR 'spgr_' num2str((ii-1), '%02.0f')];
  opts1 = [opts ' -init ' dir.COREG dir.SPGR 'spgr_' num2str((ii-1), '%02.0f') '.txt '];
  
  disp(['--Flip Angle ' num2str(alpha_spgr(ii), '%02.0f') '--']);
  % Source File
  in  = [dir.BASE dir.SPGR 'spgr_'  num2str(ii, '%02.0f')];
  % Result Filename
  out = [dir.COREG dir.SPGR 'spgr_' num2str(ii, '%02.0f')];
  
  eval(['!flirt -in ' in ' -ref ' ref ' -out ' out ' -omat ' out '.txt' opts1]);
  eval(['!fslchfiletype NIFTI ' out]);
end

% Next, set T1w SPGR as ref for IR-SPGR or AFI, since they have more similar
% contrast
ref = [dir.COREG dir.SPGR 'spgr_' num2str(length(alpha_spgr), '%02.0f')];


% Coreg The FA Calbiration Map
if afi_flag == 0
  % IR-SPGR
  disp('== Coregistration of IR-SPGR ==');
  % Source File
  in  = [dir.BASE dir.IRSPGR 'irspgr'];
  % Result Filename
  out = [dir.COREG dir.IRSPGR 'irspgr'];
  eval(['!flirt -in ' in ' -ref ' ref ' -out ' out ' -omat ' out '.txt' opts]);
  eval(['!fslchfiletype NIFTI ' out]);
else
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
  in  = [dir.BASE dir.EXT_CAL 'afi_02'];
  out = [dir.COREG dir.EXT_CAL 'afi_02'];
  eval(['!flirt -in ' in ' -ref ' ref ' -out ' out ' -applyxfm -init ' dir.COREG dir.EXT_CAL 'afi_01.txt' opts]);
  eval(['!fslchfiletype NIFTI ' out]);
end

% Apply Identity XFORM to IDEAL to Make Sure Resolution Matches
if ideal_flag == 1
  disp('== Apply I4 Transform to IDEAL B0 Map ==');
  % Write out Identity Xform
  write_xform(eye(4), [dir.COREG dir.EXT_CAL 'IDEAL.txt']);
  in  = [dir.BASE dir.EXT_CAL 'IDEAL-Omega'];
  out = [dir.COREG dir.EXT_CAL 'IDEAL-Omega'];
  eval(['!flirt -in ' in ' -ref ' ref ' -out ' out ' -applyxfm -init ' dir.COREG dir.EXT_CAL 'IDEAL.txt' opts]);
  eval(['!fslchfiletype NIFTI ' out]);
end

% Next, coregister the first SSFP-180 image to the PD-weighted SPGR
% Because they have similar contrast, then coreg all SSFP images to the
% first SSFP-180 image.

ref = [dir.COREG dir.SPGR 'spgr_01'];

disp('== Coregistration of SSFP-180 ==');
disp(['--Flip Angle ' num2str(alpha_ssfp(1), '%02.0f') ' Cycle: 180--']);
in  = [dir.BASE dir.SSFP_180 'ssfp_180_01'];
% Result Filename
out = [dir.COREG dir.SSFP_180 'ssfp_180_01'];
eval(['!flirt -in ' in ' -ref ' ref ' -out ' out ' -omat ' out '.txt' opts]);
eval(['!fslchfiletype NIFTI ' out]);

% Update Ref Image
ref = [dir.COREG dir.SSFP_180 'ssfp_180_01'];

for ii = 2:length(alpha_ssfp)
  disp(['--Flip Angle ' num2str(alpha_ssfp(ii), '%02.0f') ' Cycle: 180--']);
  in  = [dir.BASE  dir.SSFP_180 'ssfp_180_' num2str(ii, '%02.0f')];
  % Result Filename
  out = [dir.COREG dir.SSFP_180 'ssfp_180_' num2str(ii, '%02.0f')];
  eval(['!flirt -in ' in ' -ref ' ref ' -out ' out ' -omat ' out '.txt' opts]);
  eval(['!fslchfiletype NIFTI  ' out]);
end


% Then, Coreg 1st SSFP-0 image to partly T2-w SSFP-180 image (4th FA SSFP)
disp('== Coregistration of SSFP-0 Data ==');
ref = [dir.COREG dir.SSFP_180 'ssfp_180_' num2str(floor(length(alpha_ssfp)/2), '%02.0f')];
ii = 1;
disp(['--Flip Angle ' num2str(alpha_ssfp(ii), '%02.0f') ' Cycle: 0--']);
in  = [dir.BASE  dir.SSFP_0 'ssfp_0_' num2str(ii, '%02.0f')];
% Result Filename
out = [dir.COREG dir.SSFP_0 'ssfp_0_' num2str(ii, '%02.0f')];
% Add -init matrix
opts1 = [opts ' -init ' dir.COREG dir.SSFP_180 'ssfp_180_' num2str(floor(length(alpha_ssfp)/2), '%02.0f') '.txt '];
eval(['!flirt -in ' in ' -ref ' ref ' -out ' out ' -omat ' out '.txt' opts1]);
eval(['!fslchfiletype NIFTI  ' out  ]);

for ii = 2:length(alpha_ssfp)

  % Finally, coreg all the rest of the SSFP-0 images to daisy-chaing SSFP-0 Images
  ref   = [dir.COREG dir.SSFP_0 'ssfp_0_' num2str((ii-1), '%02.0f')];
  opts1 = [opts ' -init ' dir.COREG dir.SSFP_0 'ssfp_0_' num2str((ii-1), '%02.0f') '.txt '];

  disp(['--Flip Angle ' num2str(alpha_ssfp(ii), '%02.0f') ' Cycle: 0--']);
  in  = [dir.BASE  dir.SSFP_0 'ssfp_0_' num2str(ii, '%02.0f')];
  % Result Filename
  out = [dir.COREG dir.SSFP_0 'ssfp_0_' num2str(ii, '%02.0f')];

  eval(['!flirt -in ' in ' -ref ' ref ' -out ' out ' -omat ' out '.txt' opts1]);
  eval(['!fslchfiletype NIFTI  ' out]);
end

time.coreg_end = datetime();
disp(['Coreg Complete: ' time.coreg_end]);

% Coreg Flag
status.coreg = 1;

save('_mcdespot_settings', 'tr_spgr', 'tr_irspgr', 'tr_ssfp', 'flags', ...
     'alpha_spgr', 'alpha_afi', 'alpha_irspgr', 'alpha_ssfp', 'ti_irspgr', 'npe_irspgr', ...
     'info_spgr', 'info_irspgr', 'info_afi', 'info_ssfp_0', 'info_ssfp_180', 'info_ideal', 'dir', 'time', 'status');
    
diary('off');