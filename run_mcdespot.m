% FUNCTION [mcd_fv mcd_rnrm] = run_mcdespot()
%
% FUNCTION to automatie call to mcDESPOT fitting routine.
%
% Inputs:
%    EXPLICIT: None
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
% University of Oxford
% v5.5   7-Jul-2015
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
%     v5.0 - Make maj rev match between all mcDESPOT codes (Mar-2014)
%     v5.1 - Use [dir.BASE dir.SPGR 'spgr_01'] instead of info_spgr for writing out
%            new NIfTI file headers. (Avoids FSL orientation mis-label if flags.reorient
%            is set.) (Apr-2014)
%     v5.2 - Write out residual as well as maps to NIfTI
%     v5.3 - Cleanup multislice processing to skip slices without masked signal.
%            Change mcdespot_model_fit call to include number of pthreads to use
%            (Jul-2015)
%     v5.4 - Setup to resume processing if a run aborted midway. Cleanup unnecessary code (Jul-2015)

function [] = run_mcdespot()

%% Initalization
tic;

% Initiate Diary
diary('_mcdespot_log.txt');

VER = 5.4;
VERDATE = '7-Jul-2015';

% E-mail Notification When Processing is Complete
NOTIFY_EMAIL = 'shurley@fmrib.ox.ac.uk';            % Sam's FMRIB Mail

% Number of threads for parallelization of mcDESPOT residuals
numThreads = feature('NumCores');

% Builtin options
DEBUG  = 0;           % Plot data fit quality
SMOOTH = 0;           % Smooth data?

% Display banner
disp('=== cpMCDESPOT - Multicomponent Relaxomtery Analysis ===');
disp('     mcDESPOT Script       (run_mcdespot)'               );
disp('     Samuel A. Hurley      shurley@wisc.edu'             );
disp('     Pouria Mossahebi      mossahebi@wisc.edu'           );
disp(['     Version ' num2str(VER, '%01.1f') '           ' VERDATE]);
disp('     FOR USE ONLY AT UNIVERSITY OF WISCONSIN.'           );
disp('========================================================');

% Load in mcdespot settings file
load _mcdespot_settings

% Check if the processing was run before
if isfield(status, 'mcdespot') && status.mcdespot > 0
  % Do logic to finish current run vs. start new run later down.
  
else
  % First processing run - define fields
  status.mcdespot   = 1;
  status.nsliceproc = 0;
  time.mcdespot_start = [];
  time.mcdespot_end   = [];
  dir.MCDESPOT = './multiComponent/';
end

% Starting run time
time.mcdespot_start = strvcat(time.mcdespot_start, datetime_stamp()); %#ok<DSTRVCT>
disp(['Processing Run Started: ' datetime_stamp()]);

% Note: Don't save settings file unless at least 1 slice has been processed

%% Data Loading

% Try to load the T1 and flip angle map, if they exist
if isfield(status, 'despot1') && status.despot1 == 1 %#ok<*NODEF>
  % Load the t1 map
  pd_spgr = load_nifti([dir.DESPOT1 'DESPOT1-PD.nii']);
  t1      = load_nifti([dir.DESPOT1 'DESPOT1-T1.nii']);
  fam     = load_nifti([dir.DESPOT1 'DESPOT1-FAM.nii']);
  disp('DESPOT1 data loaded.');
else
  error('Must run DESPOT1-HIFI first to obtain FAM and T1 map');
end

% Try to load the T2 and off-resonance map, if they exist
if isfield(status, 'despot2') && status.despot2 == 1
  % Load the t1 map
  pd_ssfp = load_nifti([dir.DESPOT1 'DESPOT2-PD.nii']);
  t2      = load_nifti([dir.DESPOT1 'DESPOT2-T2.nii']);
  omega   = load_nifti([dir.DESPOT1 'DESPOT2-Omega.nii']);
  disp('DESPOT2 data loaded.');
else
  error('Must run DESPOT2-FM first to obtain Omega and T2 map');
end

% Load MR Data. Will automatically choose coreg data, if it exists
img = load_mcdespot_series();

spgr     = img.spgr;
ssfp_0   = img.ssfp_0;
ssfp_180 = img.ssfp_180;
clear img;

% Try to load the mask, if it exists
if isfield(status, 'mask') && status.mask == 1
  % Load the mask
  mask = load_nifti([dir.MASK status.maskname]);
  disp(['User-supplied mask loaded.']);
else
  % Threshold above 0
  mask = spgr(:,:,:,1) > 0;
  disp('Using threshold mask.');
end

% Apply mask to SPGR data
dataSize = size(spgr);
spgr     = spgr     .* repmat(mask, [1 1 1 dataSize(4)]);

% Rescale data 
data_spgr     = spgr     ./ status.despot1_signalScale;
data_ssfp_0   = ssfp_0   ./ status.despot2_signalScale;
data_ssfp_180 = ssfp_180 ./ status.despot2_signalScale;
clear spgr ssfp_0 ssfp_180;

% Smooth data, if option is selected
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

% Determine array size
sizex   = dataSize(1);
sizey   = dataSize(2);
nslices = dataSize(3);


%% Resume Processing, if Previous Run Was Incomplete

% Allocate output matrix
mcd_fv   = zeros([sizex sizey nslices 6]);
mcd_rnrm = zeros([sizex sizey nslices  ]);

% Check if previous processing run is complete
if status.nsliceproc(status.mcdespot) == 0
  % Run has not started yet. Do current run
  startSlice = 1;
  disp(['Starting mcDESPOT processing run #' num2str(status.mcdespot)]);
  
elseif status.nsliceproc(status.mcdespot) < nslices
  % Last run did not finish. Try to resume current run.
  
  % Load current data
  if exist(['mcdespot_run' num2str(status.mcdespot, '%02.0f') '.mat'], 'file');
    % Load mcd_fv and mcd_rnrm from saved file
    load(['mcdespot_run' num2str(status.mcdespot, '%02.0f') '.mat']);
    
    startSlice = status.nsliceproc(status.mcdespot) + 1;
    disp(['Resuming mcDESPOT processing run #' num2str(status.mcdespot) ' from slice #' num2str(startSlice)]);
    
  else
    % Can't find saved data, start from beginning
    status.nsliceproc(status.mcdespot) = 0;
    startSlice = 1;
    disp('Unable to load saved intermediate data file.');
    disp(['Restarting mcDESPOT processing run #' num2str(status.mcdespot) ' from beginning.']);
    
  end
  
else
  % Last run was fully finished. Start a new run from the beginning.
  status.mcdespot = status.mcdespot  + 1;  % Increment counter
  status.nsliceproc(status.mcdespot) = 0; % Make new column in slice counter
  startSlice = 1; % Start at 1st slice
  disp(['Starting mcDESPOT processing run #' num2str(status.mcdespot)]);
end


%% Loop Over Slices. Save data after each slice in case processing is aborted midway
for ii = startSlice:nslices
  
  % Grab the masked spgr slice data
  data_spgr_sl     = reshape(data_spgr(:,:,ii,:),     [sizex*sizey length(alpha_spgr)]);
  
  % If sliced is empty (mask SPGR is zero), skip
  if ~(sum(data_spgr_sl(:)) == 0)
    
    % Grab rest of slice data
    data_ssfp_0_sl   = reshape(data_ssfp_0(:,:,ii,:),   [sizex*sizey length(alpha_ssfp)]);
    data_ssfp_180_sl = reshape(data_ssfp_180(:,:,ii,:), [sizex*sizey length(alpha_ssfp)]);
    pd_spgr_sl       = pd_spgr(:,:,ii);
    pd_ssfp_sl       = pd_ssfp(:,:,ii);
    t1_sl            = t1(:,:,ii);
    t2_sl            = t2(:,:,ii);
    fam_sl           = fam(:,:,ii);
    omega_sl         = omega(:,:,ii);
    
    % Display work in progress
    disp(['Working on slice #' num2str(ii) ' / ' num2str(nslices)]);
    
    % Formulate initial guess vector
    ig = [t1_sl(:) t2_sl(:) pd_spgr_sl(:) abs(pd_ssfp_sl(:))];
    
    diary('off');
    [fv rnrm] = mcdespot_model_fit(data_spgr_sl, data_ssfp_0_sl, data_ssfp_180_sl, alpha_spgr, alpha_ssfp, tr_spgr, tr_ssfp, fam_sl(:), omega_sl(:), ig, numThreads, DEBUG);
    diary('_mcdespot_log.txt');
    
    % Reshape outputs
    mcd_fv(:,:,ii,:)   = reshape(fv,   [sizex sizey 1 6]);
    mcd_rnrm(:,:,ii)   = reshape(rnrm, [sizex sizey 1  ]);
    
    % Display the time the slice finished
    disp(['Slice Complete: ' datetime_stamp()]);
    % Increment the nsliceproc counter
    status.nsliceproc(status.mcdespot) = ii;
    
    % Save data after each slice, in case of a crash
    save(['mcdespot_run' num2str(status.mcdespot,'%02.0f')], 'mcd_fv', 'mcd_rnrm');

    save('_mcdespot_settings', 'tr_spgr', 'tr_irspgr', 'tr_ssfp', 'flags', ...
        'alpha_spgr', 'alpha_afi', 'alpha_irspgr', 'alpha_ssfp', 'ti_irspgr', 'npe_irspgr', ...
        'info_spgr', 'info_irspgr', 'info_afi', 'info_ssfp_0', 'info_ssfp_180', 'info_ideal', 'dir', 'time', 'status');
    
  end
end

%% Save final information

% End run time
time.mcdespot_end = strvcat(time.mcdespot_end, datetime_stamp()); %#ok<DSTRVCT>
disp(['Processing Run Complete: ' datetime_stamp()]);

save('_mcdespot_settings', 'tr_spgr', 'tr_irspgr', 'tr_ssfp', 'flags', ...
     'alpha_spgr', 'alpha_afi', 'alpha_irspgr', 'alpha_ssfp', 'ti_irspgr', 'npe_irspgr', ...
     'info_spgr', 'info_irspgr', 'info_afi', 'info_ssfp_0', 'info_ssfp_180', 'info_ideal', 'dir', 'time', 'status');

% Text message notification
send_mail_message(NOTIFY_EMAIL, ' run_mcdespot ', ['Processing Run For ' pwd ' Complete At: ' datetime_stamp()]);
   
%% Write out to NIfTI

% Make sure multiComponent directory exists
mkdir(dir.MCDESPOT);

% Reference image (for header info)
if status.coreg == 0
  refHdr = [dir.BASE  dir.SPGR 'spgr_01.nii'];
elseif status.coreg == 1
  refHdr = [dir.COREG dir.SPGR 'spgr_01.nii'];
end

% Save NIfTI
img_nifti_to_nifti(mcd_fv(:,:,:,1), refHdr, [dir.MCDESPOT 'mcDESPOT-T1m' ]);
img_nifti_to_nifti(mcd_fv(:,:,:,2), refHdr, [dir.MCDESPOT 'mcDESPOT-T1f' ]);
img_nifti_to_nifti(mcd_fv(:,:,:,3), refHdr, [dir.MCDESPOT 'mcDESPOT-T2m' ]);
img_nifti_to_nifti(mcd_fv(:,:,:,4), refHdr, [dir.MCDESPOT 'mcDESPOT-T2f' ]);
img_nifti_to_nifti(mcd_fv(:,:,:,5), refHdr, [dir.MCDESPOT 'mcDESPOT-MWF' ]);
img_nifti_to_nifti(mcd_fv(:,:,:,6), refHdr, [dir.MCDESPOT 'mcDESPOT-Tau' ]);
img_nifti_to_nifti(mcd_rnrm       , refHdr, [dir.MCDESPOT 'mcDESPOT-Rnrm']);
     
diary('off');
