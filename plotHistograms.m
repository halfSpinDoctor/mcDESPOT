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
figure;

% Make Subplots
subplot(2,3,1);
hist(mwf(find(mwf)),100);
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
figure;

% Plot MWF parameters
subplot(2,3,1);
imagesc(mwf(:,:,sl), [0 .30]);
axis image; axis off; colormap gray; colorbar;
title 'MWF'


subplot(2,3,2)
imagesc(t1m(:,:,sl));
axis image; axis off; colormap gray; colorbar;
title 'T1_m'
xlabel 's'

subplot(2,3,3)
imagesc(t1f(:,:,sl));
axis image; axis off; colormap gray; colorbar;
title 'T1_f'
xlabel 's'

subplot(2,3,4)
imagesc(tau(:,:,sl));
axis image; axis off; colormap gray; colorbar;
title 'Tau'
xlabel 's'

subplot(2,3,5)
imagesc(t2m(:,:,sl));
axis image; axis off; colormap gray; colorbar;
title 'T2_m'
xlabel 's'

subplot(2,3,6)
imagesc(t2f(:,:,sl));
axis image; axis off; colormap gray; colorbar;
title 'T2_f'
xlabel 's'