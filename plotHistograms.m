%% FUNCTION [] = plotHistograms()
%
% Function to plot images and histograms of the multicomponent parameter maps in the
% current directory (assumes mcDESPOT- prefix and NIfTI files)
%
% Samuel A. Hurley
% University of Wisconsin
% v1.0 5-Feb-2014
%
% Changelog:
%   v1.0 - Initial version (Feb 2014)

% Load Data
mwf = load_nifti('mcDESPOT-MWF.nii');
mwf = load_nifti('mcDESPOT-MWF.nii');
t1f = load_nifti('mcDESPOT-T1f.nii');
t1m = load_nifti('mcDESPOT-T1m.nii');
t2f = load_nifti('mcDESPOT-T2f.nii');
t2m = load_nifti('mcDESPOT-T2m.nii');
tau = load_nifti('mcDESPOT-Tau.nii');

% Trim near-zero values for MWF
mwf(mwf<.005) = 0;

% Make Subplots
subplot(2,3,1);
title 'MWF'
xlabel '%'
subplot(2,3,2)
hist(t1m(find(t1m)),100);
title 'T1_m'
xlabel 's'
subplot(2,3,3)
hist(t1f(find(t1f)),100);
title 'T1_f'
xlabel 's'
subplot(2,3,4)
hist(tau(find(tau)),100);
title 'Tau'
xlabel 's'
subplot(2,3,5)
hist(t2m(find(t2m)),100);
title 'T2_m'
xlabel 's'
subplot(2,3,6)
hist(t2f(find(t2f)),100);
title 'T2_f'
xlabel 's'

% Slice for plotting
sl = 84;

% Plot MWF parameters
subplot(2,3,1);
imagesc(mcd_fv(:,:,sl,5));
axis image; axis off;
title 'MWF'
xlabel '%'

subplot(2,3,2)
imagesc(mcd_fv(:,:,sl,1));
axis image; axis off;
title 'T1_m'
xlabel 's'

subplot(2,3,3)
imagesc(mcd_fv(:,:,sl,2));
axis image; axis off;
title 'T1_f'
xlabel 's'

subplot(2,3,4)
imagesc(mcd_fv(:,:,sl,6));
axis image; axis off;
title 'Tau'
xlabel 's'

subplot(2,3,5)
imagesc(mcd_fv(:,:,sl,3));
axis image; axis off;
title 'T2_m'
xlabel 's'

subplot(2,3,6)
imagesc(mcd_fv(:,:,sl,4));
axis image; axis off;
title 'T2_f'
xlabel 's'