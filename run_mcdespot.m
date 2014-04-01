% FUNCTION [mcd_fv mcd_rnrm] = run_mcdespot(slices)
%
% FUNCTION to automatie call to mcDESPOT fitting routine.
%
% Inputs:
%    EXPLICIT: slices - choose which slices to process
%    IMPLICIT: Start in same directory as _mcdespot_settings.mat
%              Have already run brain mask, T1, and FAM steps (despot1)
%
% Outputs:
%    EXPLICIT: None
%    IMPLICIT: New directory called multiComponent
%
% Outputs:
%    mcd_fv   - "function values" 4D matrix with [t1m t1f t2m t2f Fm Tau Omega]
%                               along 4th dimension
%
%    mcd_rnrm - 3d matrix of residuals
%
% Samuel A. Hurley
% University of Wisconsin
% v4.0 09-Feb-2011
%
% Changelog:
%     v1.0 - added options for coregistered data               (Jan-2010)
%          - now saves data automatically to _mcdespot_proc_mode#
%          - now saves output data in same matrix as input (even if less
%            slices are being processed)
%          - deletes slice data from memory as they are no longer needed
%     v2.0 - Many new cool options. Calls v2.0 of fitting routine  (26-Feb-2010)
%     v3.0 - Updated to work with newest mcdespot_model_fit version 3.0 (Jun-2010)
%     v4.0 - Updated to work with newest mcdespot_model_fit version 4.2  (Feb-2011)
%            Removed all inputs and outputs to the function
%     v4.1 - Runs all slices if slices is not specified

function [] = run_mcdespot(slices)
tic;

% Initiate Diary
diary('_mcdespot_log.txt');
clc;

VER = 4.0;
VERDATE = '09-Feb-2011';

% E-mail Notification When Processing is Complete
NOTIFY_EMAIL = '2628943063@email.uscc.net';            % Sam's Cell Phone

% Builtin options
DEBUG  = 0;           % Plot data fit quality
SMOOTH = 0;           % Smooth data?

% Display banner
disp('=== cpMCDESPOT - Multicomponent Relaxomtery Analysis ===');
disp('     mcDESPOT Script       (run_mcdespot)'               );
disp('     Samuel A. Hurley      shurley@wisc.edu'             );
disp('     Pouria Mossahebi      mossahebi@wisc.edu'           );
disp(['     Version ' num2str(VER, '%01.1f') '         ' VERDATE]);
disp('     FOR USE ONLY AT UNIVERSITY OF WISCONSIN.'           );
disp('========================================================');

% Load in mcdespot settings file
load _mcdespot_settings


% Check if the processing was run before
if isfield(status, 'mcdespot') && status.mcdespot > 0 %#ok<NODEF>
  status.mcdespot        = status.mcdespot + 1;
else
  status.mcdespot = 1;
  status.comment  = [];
  time.mcdespot_start = [];
  time.mcdespot_end   = [];
  
  % Define multiComponant location
  dir.MCDESPOT = './multiComponent/';
  mkdir(dir.MCDESPOT);
end

disp(['Processing Run #' num2str(status.mcdespot)]);

% Comment for this run
% comment = inputdlg('Enter a comment for this run:');
comment = ' ';

% Starting time for processing
time.mcdespot_start = strvcat(time.mcdespot_start, datetime());  %#ok<REMFF1>
disp(['Processing Run Started: ' datetime()]);

% Load mcDESPOT Images
% load_mcdespot_series will choose the coreg data, if it exists
img = load_mcdespot_series();

% Grab the necessary data
spgr     = img.spgr;
ssfp_0   = img.ssfp_0;
ssfp_180 = img.ssfp_180;

% Free up memory
img = []; %#ok<NASGU>

% Try to load the T1 and flip angle map, if they exist
if isfield(status, 'despot1') && status.despot1 == 1
  % Load the t1 map
  pd_spgr = load_nifti([dir.DESPOT1 'DESPOT1-PD.nii']);
  t1      = load_nifti([dir.DESPOT1 'DESPOT1-T1.nii']);
  fam     = load_nifti([dir.DESPOT1 'DESPOT1-FAM.nii']);
else
  error('Must run DESPOT1-HIFI first to obtain FAM and T1 map');
end

% Try to load the T2 and off-resonance map, if they exist
if isfield(status, 'despot2') && status.despot2 == 1
  % Load the t1 map
  pd_ssfp = load_nifti([dir.DESPOT1 'DESPOT2-PD.nii']);
  t2      = load_nifti([dir.DESPOT1 'DESPOT2-T2.nii']);
  omega   = load_nifti([dir.DESPOT1 'DESPOT2-Omega.nii']);
else
  error('Must run DESPOT2-FM first to obtain Omega and T2 map');
end

% Try to load the mask, if it exists
if isfield(status, 'mask') && status.mask == 1
  % Load the mask
  mask = load_nifti([dir.MASK status.maskname]);
  disp(['Using user-supplied mask.']);
else
  % Threshold above 0
  mask = spgr(:,:,:,1) > 0;
  disp('Using threshold mask.');
end

% Mask data (SPGR Only)
dataSize = size(spgr);
spgr     = spgr     .* repmat(mask, [1 1 1 dataSize(4)]);

% Normalize data 
spgr     = spgr     ./ status.despot1_signalScale;
ssfp_0   = ssfp_0   ./ status.despot2_signalScale;
ssfp_180 = ssfp_180 ./ status.despot2_signalScale;

% Check if slices field was specified
if ~check_var('slices')
  slices = 1:size(spgr,3);
  
end

% Select slices for all maps
data_spgr = spgr(:,:,slices,:);
data_ssfp_0 = ssfp_0(:,:,slices,:);
data_ssfp_180 = ssfp_180(:,:,slices,:);

pd_spgr = pd_spgr(:,:,slices);
pd_ssfp = pd_ssfp(:,:,slices);
t1      = t1(:,:,slices);
t2      = t2(:,:,slices);
fam     = fam(:,:,slices);
omega   = omega(:,:,slices);

% Clear unused vars
clear spgr ssfp_0 ssfp_180;

% Smooth data
if SMOOTH == 1
  fprintf('Smooth Data /w sigma = 2 gaussian (5x5 kernel)...');
  fprintf('SPGR...');
  data_spgr = img_smooth(data_spgr, 5, 2);
  fprintf('SSFP-0...');
  data_ssfp_0 = img_smooth(data_ssfp_0, 5, 2);
  fprintf('SSFP-180...');
  data_ssfp_180 = img_smooth(data_ssfp_180, 5, 2);
  disp('done.');
end

% Display masked slices that will be used for processing
% imgsc(data_spgr(:,:,:,end));
% title 'Images to Process';

%% Save the run status before processing begins
status.comment                     = 'Standard mcDESPOT Processing'; %#ok<VCAT>
status.nsliceproc(status.mcdespot) = 0;

% Determine array size
sizex   = size(data_spgr, 1);
sizey   = size(data_spgr, 2);
nslices = size(data_spgr, 3);

% Create Output Matrix
mcd_fv   = zeros([sizex sizey nslices 6]);
mcd_rnrm = zeros([sizex sizey nslices  ]);

%% Loop Over Slices, so if it crashes we don't loose all the data
for ii = 1:nslices
  disp(['Working on slice #' num2str(ii) ' / ' num2str(nslices)]);
  
  % Grab top slice off of stack
  data_spgr_sl     = reshape(data_spgr(:,:,1,:),     [sizex*sizey length(alpha_spgr)]);
  data_ssfp_0_sl   = reshape(data_ssfp_0(:,:,1,:),   [sizex*sizey length(alpha_ssfp)]);
  data_ssfp_180_sl = reshape(data_ssfp_180(:,:,1,:), [sizex*sizey length(alpha_ssfp)]);
  pd_spgr_sl       = pd_spgr(:,:,1);
  pd_ssfp_sl       = pd_ssfp(:,:,1);
  t1_sl            = t1(:,:,1);
  t2_sl            = t2(:,:,1);
  fam_sl           = fam(:,:,1);
  omega_sl         = omega(:,:,1);
  
  % Remove top slice to free up memory, then the next slice will be processed next
  data_spgr(:,:,1,:)     = [];
  data_ssfp_0(:,:,1,:)   = [];
  data_ssfp_180(:,:,1,:) = [];
  pd_spgr(:,:,1)         = [];
  pd_ssfp(:,:,1)         = [];
  t1(:,:,1)              = [];
  t2(:,:,1)              = [];
  fam(:,:,1)             = [];
  omega(:,:,1)           = [];
  
  % Formulate initial guess vector
  ig = [t1_sl(:) t2_sl(:) pd_spgr_sl(:) abs(pd_ssfp_sl(:))];
  
  diary('off');
  [fv rnrm] = mcdespot_model_fit(data_spgr_sl, data_ssfp_0_sl, data_ssfp_180_sl, alpha_spgr, alpha_ssfp, tr_spgr, tr_ssfp, fam_sl(:), omega_sl(:), ig, DEBUG);
  diary('_mcdespot_log.txt');
  
  % Reshape outputs
  slicenum = ii + (slices(1) - 1);
  mcd_fv(:,:,slicenum,:)   = reshape(fv,   [sizex sizey 1 6]);      
  mcd_rnrm(:,:,slicenum)   = reshape(rnrm, [sizex sizey 1  ]);        
  
  % Display the time the slice finished
  disp(['Slice Complete: ' datetime()]);
  
  % Save data after each slice, in case of a crash
  save(['mcdespot_run' num2str(status.mcdespot,'%02.0f')], 'mcd_fv', 'mcd_rnrm');
  status.nsliceproc(status.mcdespot) = status.nsliceproc(status.mcdespot) + 1;
  
  save('_mcdespot_settings', 'tr_spgr', 'tr_irspgr', 'tr_ssfp', 'flags', ...
     'alpha_spgr', 'alpha_afi', 'alpha_irspgr', 'alpha_ssfp', 'ti_irspgr', 'npe_irspgr', ...
     'info_spgr', 'info_irspgr', 'info_afi', 'info_ssfp_0', 'info_ssfp_180', 'info_ideal', 'dir', 'time', 'status');
    
end

% Save end run time
time.mcdespot_end = strvcat(time.mcdespot_end, datetime()); %#ok<VCAT>
disp(['Processing Run Complete: ' datetime()]);

% Text message notification
send_mail_message(NOTIFY_EMAIL, ' run_mcdespot ', ['Processing Run Complete: ' datetime()]);

%% Write out to NIfTI
% Save NIfTI
img_dcm_to_nifti(mcd_fv(:,:,:,1), info_spgr,[dir.MCDESPOT 'mcDESPOT-T1m']);
img_dcm_to_nifti(mcd_fv(:,:,:,2), info_spgr,[dir.MCDESPOT 'mcDESPOT-T1f']);
img_dcm_to_nifti(mcd_fv(:,:,:,3), info_spgr,[dir.MCDESPOT 'mcDESPOT-T2m']);
img_dcm_to_nifti(mcd_fv(:,:,:,4), info_spgr,[dir.MCDESPOT 'mcDESPOT-T2f']);
img_dcm_to_nifti(mcd_fv(:,:,:,5), info_spgr,[dir.MCDESPOT 'mcDESPOT-MWF']);
img_dcm_to_nifti(mcd_fv(:,:,:,6), info_spgr,[dir.MCDESPOT 'mcDESPOT-Tau']);


%%
save('_mcdespot_settings', 'tr_spgr', 'tr_irspgr', 'tr_ssfp', 'flags', ...
     'alpha_spgr', 'alpha_afi', 'alpha_irspgr', 'alpha_ssfp', 'ti_irspgr', 'npe_irspgr', ...
     'info_spgr', 'info_irspgr', 'info_afi', 'info_ssfp_0', 'info_ssfp_180', 'info_ideal', 'dir', 'time', 'status');
     
diary('off');