% FUNCTION [] = normalize_mcdespot();
%
% FUNCTION to apply spatial warping transform computed from ANTs to individual
%          mcDESPOT datasets for population-level / voxelwise analysis
%
% Inputs:
%    Explicit:
%              NONE
%
%    Implicit: Start in same directory as _mcdespot_settings.mat
%
% Outputs:
%    Explicit: None
%    Implicit: New directory called normalizedData
%
% Samuel A. Hurley
% v5.0 - 22-Jul-2014
%
% Changelog:
%     v5.0 - Initial version (start at 5 to match major version of other funcs)

function normalize_mcdespot()

% Diary
diary('_mcdespot_log.txt');

% Load mcDESPOT Settings
load '_mcdespot_settings';

% Time
time.norm_start = datetime_stamp;
disp(['Spatial Normalization Started: ' time.norm_start]);

% Define Spatial Normalization Directory
dir.NORM   = './normalizedData/';
mkdir(dir.NORM);
mkdir([dir.NORM dir.DESPOT1]);
mkdir([dir.NORM dir.MCDESPOT]);

% Define location of ANTs binaries
ANTSBIN    = '/home/samuel/Software/ANTs/bin/';

% Define Transforms To Apply
ants_xform = {'ANTs_Warp.nii.gz' 'ANTs_Affine.txt'};

seriesOfXforms = [];
for jj = 1:length(ants_xform)
  seriesOfXforms = horzcat(seriesOfXforms, [dir.NORM ants_xform{jj} ' ']); %#ok<AGROW>
end

% Apply Transforms to Maps

% DESPOT1
if isfield(status, 'despot1') && status.despot1 == 1 %#ok<NODEF>
  maps = {'DESPOT1-FAM.nii' 'DESPOT1-FAM_Unsmooth.nii' 'DESPOT1-PD.nii' 'DESPOT1-Rnrm.nii' 'DESPOT1-T1.nii'};
  
  for ii = 1:length(maps)
    in    = [dir.DESPOT1 maps{ii}];
    out   = [dir.NORM    dir.DESPOT1 maps{ii}];

    cmd   = ['!' ANTSBIN 'WarpImageMultiTransform 3 ' in ' ' out ' ' seriesOfXforms];
    eval(cmd);
  end
end

% DESPOT2
if isfield(status, 'despot2') && status.despot2 == 1
  maps = {'DESPOT2-Omega.nii' 'DESPOT2-PD.nii' 'DESPOT2-Rnrm.nii' 'DESPOT2-T2.nii'};
  
  for ii = 1:length(maps)
    in    = [dir.DESPOT1 maps{ii}];
    out   = [dir.NORM    dir.DESPOT1 maps{ii}];

    cmd   = ['!' ANTSBIN 'WarpImageMultiTransform 3 ' in ' ' out ' ' seriesOfXforms];
    eval(cmd);
  end
end

% mcDESPOT
if isfield(status, 'mcdespot') && status.mcdespot == 1
  maps = {'mcDESPOT-MWF.nii' 'mcDESPOT-T1f.nii' 'mcDESPOT-T1m.nii' 'mcDESPOT-T2f.nii' 'mcDESPOT-T2m.nii' 'mcDESPOT-Tau.nii'};
  
  for ii = 1:length(maps)
    in    = [dir.MCDESPOT maps{ii}];
    out   = [dir.NORM     dir.MCDESPOT maps{ii}];

    cmd   = ['!' ANTSBIN 'WarpImageMultiTransform 3 ' in ' ' out ' ' seriesOfXforms];
    eval(cmd);
  end
end

% Complete!
status.norm = 1;
time.norm_end = datetime_stamp;
disp(['Spatial Normalization Finished: ' time.norm_end]);

save('_mcdespot_settings', 'tr_spgr', 'tr_irspgr', 'tr_ssfp', 'flags', ...
     'alpha_spgr', 'alpha_afi', 'alpha_irspgr', 'alpha_ssfp', 'ti_irspgr', 'npe_irspgr', ...
     'info_spgr', 'info_irspgr', 'info_afi', 'info_ssfp_0', 'info_ssfp_180', 'info_ideal', 'dir', 'time', 'status');
    
diary('off');