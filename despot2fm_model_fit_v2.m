% Function [pd r2 omega rnrm] = despot2fm_model_fit_v2(data_0, data_180, alpha, tr, r1, pd, fam, opts)
%
% Function to fit single-componant DESPOT2 model, including off-resonance effects
% and extra term to account for TE
%
% NOTE: Only requires data_0 to be masked
%
% Inputs:
%   data_0    - (Num Points) x (Num of Flip Angles) SSFP /w   0 phase cycling
%   data_180  - (Num Points) x (Num of Flip Angles) SSFP /w 180 phase cycling
%   alpha     - a vector of flip angles [degrees]
%   tr        - repitition time [s]
%   r1        - Spin-lattice relaxation map, as fit from DESPOT1
%   pd        - Proton density relaxation map, as fit from DESPOT1, for initial guess
%   fam       - flip angle map to account for B1 inhomogeniety
%
%   opts  - a structure containing optional settings
%       points   - use a subset of points for 
%       debug    - 0 = debug by plotting data pts and fit curve
%                  1 = debug by showing GA fitness plots & plot of curve
%
% Outputs:
%   pd    - Proton density [a.u.]
%   r2    - Spin-spin relaxation rate [s^-1]
%   omega - B0 off-resonance [Hz]
%   rnrm  - Sum-of-squares residual
%
% References:
%   DESPOT2-FM: 
%
% Samuel A. Hurley
% University of Wisconsin
% v5.2 8-Jul-2015
%
% Chagelog:
%      v1.0 - Initial version based on despot1hifi_model_fit     (Feb-2011)
%      v2.0 - Switched to GA (Genetic Algorithm) fitting routine (Jun-2011)
%      v2.1 - Normalize by mean of MR signal to avoid explicit fitting for PD in GA, similar to Sean's mcDESPOT (Jun-2011)
%      v2.2 - Two-pass process to fit first for Omega ig (assuming fixed T2), 
%             then again for T2, PD, and Omega. Option for external B0 map  (Oct-2011)
%      v2.3 - Input PD initial guess from DESPOT1. Use R1 to formulate T2
%             initial guess based on Vasily constraints. (Nov-2011)
%      v2.4 - Removed unused GA stuff, updated pd initial guess for ext. supplied B0 map (Jan-2012)
%      v5.0 - Update to be compatible /w new cpDESPOT2 C-Code (Nov-2012)
%      v5.1 - Implement version with parfor parallelisation, similar to
%             despot1hifi_model_fit_v2 (Jun-2015)
%      v5.2 - Embed median filtering function (medfilt3) into this m-file (Jul-2015)

function [pd r2 omega rnrm] = despot2fm_model_fit_v2(data_0, data_180, alpha, tr, r1, pd_despot1, fam, opts)

tic;

% Constant: Assumption of ratio of T2 to T1 for 1st pass
R1_R2_RATIO = 0.045; % Yarhykh assumption is 0.045

% Nelder-Mead Downhill Simplex (fminsearch)
optim=optimset('fminsearch');
optim.TolFun = 1e-4;
optim.TolX   = 1e-4;
optim.MaxFunEvals = 100;
optim.Display = 'off';

% I.Check input arguments
switch nargin
  case 7
    % No options struct specified, use default options
    opts   = struct();    % Use default options
  case 8
    if ~isstruct(opts)
      error('Opts must be a structure.  Type ''help despot2fm_model_fit'' for more information');
    end
  otherwise
    error('You must supply 7 or 8 input arguments.  Type ''help despot2fm_model_fit'' for more information');
end

% Grab Options or Set Defaults
if isfield(opts, 'points')
  points = opts.points;
end

if isfield(opts, 'debug')
  debug = opts.debug;
  if debug > 0
    dbgtime = toc;
  end
else
  % Default is to turn off debugging
  debug = 0;
end

if isfield(opts, 'ext_omega')
  % Externally supplied B0 map, in Hz. -- NOT SUPPORTED FOR PARFOR VERSION
  error('External B0 map is not supported for parallel version');
end


% Banner
disp('==============DESPOT2-FM T2===================');

% II. Check data dimensions
if size(alpha, 2) ~= size(data_0, 2)
  error('Number of SSFP-0 flip angles and supplied SPGR data points do not match.');
end

if size(alpha, 2) ~= size(data_180, 2)
  error('Number of SSFP-180 flip angles and supplied SPGR data points do not match.');
end

% Determine which points to use
if exist('points', 'var')
  % Use specified points
  if length(points) > length(alpha)
    error('You specifed more fitting points than avalible flip angles');
  end
  
  alpha = alpha(points);
  disp('NOTE: Subset of SSFP flip angles used for fitting.');
  
  % Grab only the flip angles specified
  data_0   = data_0(:,points);
  data_180 = data_180(:,points);
end

disp(['SSFP Flip Angles: ' num2str(alpha, '%01.0f ')]);

% Apply Flip Angle Correction
alpha = fam * alpha;

% Convert to radians
alpha    = alpha .* pi/180;

% Preallocate output [pd r2 omega]
pd    = zeros([size(data_0, 1) 1]);
r2    = zeros([size(data_0, 1) 1]);
omega = zeros([size(data_0, 1) 1]);
rnrm  = zeros([size(data_0, 1) 1]);

npts = size(find(~(sum(data_0, 2) == 0)),1);

% Preallocate parfor output [pd r2 omega]
pd_pf    = zeros([npts 1]);
r2_pf    = zeros([npts 1]);
omega_pf = zeros([npts 1]);
rnrm_pf  = zeros([npts 1]);

fprintf('DESPOT2-FM Fitting:...');

% III. Perform DESPOT2-FM Fit in Each Voxel

% Get non-zero voxels
voxidx = find(~(sum(data_0, 2) == 0))';

parfor ii = 1:length(voxidx)
  
  % Get parfor index
  ii_pf = voxidx(ii);
  
  % Grab voxel data
  vox_data_0   = data_0(ii_pf,:);     %#ok<*PFBNS>
  vox_data_180 = data_180(ii_pf,:);
  
  vox_alpha    = alpha(ii_pf, :);
  vox_r1       = r1(ii_pf);
  vox_pd_d1    = pd_despot1(ii_pf);
  
%   % Update Progress
%   <<REOMVED>>
  
  % Minimize with Fminsearch
  
  % -- FIRST PASS: fix T2, fit for initial guess of Omega --
  % Fit for PD and Off-Res: fix_flag = 1
  ig = [vox_pd_d1 10];                     % PD / Omega initial guess
  r2_ig = vox_r1 / R1_R2_RATIO;            % Fix T2 based on T2/T1 ~= 0.045 (0.055 from Yarnykh NI2004)
  despot2fm_model_handle_1 = @(x)(cpDESPOT2_residuals_SAH([x(1) 1./vox_r1 1./r2_ig x(2)], 0, vox_data_0', vox_alpha', tr, 1)' + cpDESPOT2_residuals_SAH([x(1) 1./vox_r1 1./r2_ig x(2)], 180, vox_data_180', vox_alpha', tr, 1)');
  [x1] = fminsearch(despot2fm_model_handle_1, ig, optim);

  % -- SECOND PASS - set omega IG from previous step --
  ig = [x1(1) vox_r1/R1_R2_RATIO x1(2)];   % PD from X1, R2 = R1/0.045, Omega from X1
  despot2fm_model_handle_2 = @(x)(cpDESPOT2_residuals_SAH([x(1) 1./vox_r1 1./x(2)  x(3)], 0, vox_data_0', vox_alpha', tr, 1)' + cpDESPOT2_residuals_SAH([x(1) 1./vox_r1 1./x(2)  x(3)], 180, vox_data_180', vox_alpha', tr, 1)');
  [x2]   = fminsearch(despot2fm_model_handle_2, ig, optim);
  
  % Save results
  pd_pf(ii)    = x2(1);
  r2_pf(ii)    = x2(2);
  omega_pf(ii) = x2(3);
  
end

% Expand output variables
pd(voxidx)    = pd_pf;
r2(voxidx)    = r2_pf;
omega(voxidx) = omega_pf;

fprintf('First 2 iterations complete...applying median filter...');

% Apply 3x3x3 3D median filter to Omega estimate
omega_mf5     = medfilt3(omega, [7 7 7]);

fprintf('median filtering complete...running last iteration...');


parfor ii = 1:length(voxidx)
  
  % Get parfor index
  ii_pf = voxidx(ii);
  
  % Grab voxel data
  vox_data_0   = data_0(ii_pf,:);
  vox_data_180 = data_180(ii_pf,:);
  
  vox_alpha    = alpha(ii_pf, :);
  vox_r1       = r1(ii_pf);
  
  vox_pd_ig    = pd(ii_pf);
  vox_r2_ig    = r2(ii_pf);
  vox_omega_ig = omega_mf5(ii_pf);

  % -- THIRD PASS - re-run second pass, with median filtered omega as ig --
  ig = [vox_pd_ig vox_r2_ig vox_omega_ig];    % Initial guess: PD from X2, R2 from X2, Omega from median filtered image
  despot2fm_model_handle_3 = @(x)(cpDESPOT2_residuals_SAH([x(1) 1./vox_r1 1./x(2) x(3)], 0, vox_data_0', vox_alpha', tr, 1)' + cpDESPOT2_residuals_SAH([x(1) 1./vox_r1 1./x(2) x(3)], 180, vox_data_180', vox_alpha', tr, 1)');
  [x residual]   = fminsearch(despot2fm_model_handle_3, ig, optim);
  
  % Save results
  pd_pf(ii)    = x(1);
  r2_pf(ii)    = x(2);
  omega_pf(ii) = x(3);
  rnrm_pf(ii)  = residual;
  
%   % PLOT FOR DEBUGGING FINAL FIT RESULT
%   <<REMOVED>>
  
end

% Close parfor progress bar
% parfor_progress(0);

% Expand output variables
pd(voxidx)    = pd_pf;
r2(voxidx)    = r2_pf;
omega(voxidx) = omega_pf;
rnrm(voxidx)  = rnrm_pf;

% Done.
disp('...fitting complete.');
toc;

%% == Helper Functions Below ==

function B = medfilt3(A,siz,padopt,CHUNKFACTOR)

%MEDFILT3 1-D, 2-D and 3-D median filtering.
%   B = MEDFILT3(A,[M N P]) performs median filtering of the 3-D array A.
%   Each output pixel contains the median value in the M-by-N-by-P
%   neighborhood around the corresponding pixel in the input array.
%
%   B = MEDFILT3(A,[M N]) performs median filtering of the matrix A. Each
%   output pixel contains the median value in the M-by-N neighborhood
%   around the corresponding pixel.
%
%   B = MEDFILT3(A,M) performs median filtering of the vector A. Each
%   output pixel contains the median value in the M neighborhood
%   around the corresponding pixel.
%
%   B = MEDFILT3(A) performs median filtering using a 3 or 3x3 or 3x3x3
%   neighborhood according to the size of A.
%
%   B = MEDFILT3(A,...,PADOPT) pads array A using PADOPT option:
%
%      String values for PADOPT (default = 'replicate'):
%      'circular'    Pads with circular repetition of elements.
%      'replicate'   Repeats border elements of A. (DEFAULT)
%      'symmetric'   Pads array with mirror reflections of itself.
%
%      If PADOPT is a scalar, A is padded with this scalar.
%
%   Class Support
%   -------------
%     Input array can be numeric or logical. The returned array is of class
%     single or double.
%
%   Notes
%   -----
%     M, N and P must be odd integers. If not, they are incremented by 1.
%
%     If NANMEDIAN exists (Statistics Toolbox is required), then MEDFILT3
%     treats NaNs as missing values.
%
%     If you work with very large 3D arrays, an "Out of memory" error may
%     appear. The chunk factor (CHUNKFACTOR, default value = 1) must be
%     increased to reduce the size of the chunks. This will imply more
%     iterations whose number is directly proportional to CHUNKFACTOR. Use
%     the following syntax: MEDFILT3(A,[...],PADOPT,CHUNKFACTOR)
%
%   Examples
%   --------
%     %>> 1-D median filtering <<
%     t = linspace(0,2*pi,100);
%     y = cos(t);
%     I = round(rand(1,5)*99+1);
%     y(I) = rand(size(I));
%     ys = medfilt3(y,5);
%     plot(t,y,':',t,ys)
%
%     %>> 2-D median filtering <<
%     % original image
%     I = imread('eight.tif');
%     % noisy image
%     J = I;
%     rand('state',sum(100*clock))
%     J(rand(size(J))<0.01) = 255;
%     J(rand(size(J))<0.01) = 0;
%     % denoised image
%     K = medfilt3(J);
%     % figures
%     figure
%     subplot(121),imshow(J), subplot(122), imshow(K)
%
%     %>> 3-D median filtering <<
%     rand('state',0)
%     [x,y,z,V] = flow(50);
%     noisyV = V + 0.1*double(rand(size(V))>0.95);
%     clear V
%     figure
%     subplot(121)
%     hpatch = patch(isosurface(x,y,z,noisyV,0));
%     isonormals(x,y,z,noisyV,hpatch)
%     set(hpatch,'FaceColor','red','EdgeColor','none')
%     daspect([1,4,4]), view([-65,20]), axis tight off
%     camlight left; lighting phong 
%     subplot(122)
%     %--------
%     denoisedV = medfilt3(noisyV,7);
%     %--------
%     hpatch = patch(isosurface(x,y,z,denoisedV,0));
%     isonormals(x,y,z,denoisedV,hpatch)
%     set(hpatch,'FaceColor','red','EdgeColor','none')
%     daspect([1,4,4]), view([-65,20]), axis tight off
%     camlight left; lighting phong
%       
%   See also MEDFILT1, MEDFILT2, HMF.
%
%   -- Damien Garcia -- 2007/08, revised 2008/10

%% Note:
% If you work with large 3D arrays, an "Out of memory" error may appear.
% The chunk factor thus must be increased to reduce the size of the chunks.
if nargin~=4
    CHUNKFACTOR = 1;
end
if CHUNKFACTOR<1, CHUNKFACTOR = 1; end

%% Checking input arguments
if isscalar(A), B = A; return, end

if ndims(A)>3
    error('A must be a 1-D, 2-D or 3-D array.')
end

if all(isnan(A(:))), B = A; return, end

sizA = size(A);
if nargin==1
    % default kernel size is 3 or 3x3 or 3x3x3
    if isvector(A)
        siz = 3;
    else
        siz = 3*ones(1,numel(sizA));
    end
    padopt = 'replicate';
elseif nargin==2
    % default padding option is "replicate"
    padopt = 'replicate';
end

%% Make SIZ a 3-element array
if numel(siz)==2
    siz = [siz 1];
elseif isscalar(siz)
    if sizA(1)==1
        siz = [1 siz 1];
    else
        siz = [siz 1 1];
    end
end

%% Chunks: the numerical process is split up in order to avoid large arrays
N = numel(A);
siz = ceil((siz-1)/2);
n = prod(siz*2+1);
if n==1, B = A; return, end
nchunk = (1:ceil(N/n/CHUNKFACTOR):N);
if nchunk(end)~=N, nchunk = [nchunk N]; end

%% Change to double if needed
class0 = class(A);
if ~isa(A,'float')
    A = double(A);
end

%% Padding along specified direction
% If PADARRAY exists (Image Processing Toolbox), this function is used.
% Otherwise the array is padded with scalars.
B = A;
sizB = sizA;
try
    A = padarray(A,siz,padopt);
catch
    if ~isscalar(padopt)
        padopt = 0;
        warning('MATLAB:medfilt3:InexistentPadarrayFunction',...
            ['PADARRAY function does not exist: '...
            'only scalar padding option is available.\n'...
            'If not specified, the scalar 0 is used as default.']);
    end
    A = ones(sizB+siz(1:ndims(B))*2)*padopt;
    A(siz(1)+1:end-siz(1),siz(2)+1:end-siz(2),siz(3)+1:end-siz(3)) = B;
end
sizA = size(A);

if numel(sizB)==2
    sizA = [sizA 1];
    sizB = [sizB 1];
end

%% Creating the index arrays (INT32)
inc = zeros([3 2*siz+1],'int32');
siz = int32(siz);
[inc(1,:,:,:) inc(2,:,:,:) inc(3,:,:,:)] = ndgrid(...
    [0:-1:-siz(1) 1:siz(1)],...
    [0:-1:-siz(2) 1:siz(2)],...
    [0:-1:-siz(3) 1:siz(3)]);
inc = reshape(inc,[1 3 prod(2*single(siz)+1)]);

I = zeros([sizB 3],'int32');
sizB = int32(sizB);
[I(:,:,:,1) I(:,:,:,2) I(:,:,:,3)] = ndgrid(...
    (1:sizB(1))+siz(1),...
    (1:sizB(2))+siz(2),...
    (1:sizB(3))+siz(3));
I = reshape(I,[prod(single(sizB)) 3]);

%% Check if NANMEDIAN exists
existNaNmedian = exist('nanmedian','file');

%% Filtering
for i = 1:length(nchunk)-1

    Im = repmat(I(nchunk(i):nchunk(i+1),:),[1 1 n]);
    Im = Im + repmat(inc,[nchunk(i+1)-nchunk(i)+1,1,1]);

    I0 = Im(:,1,:) +...
        (Im(:,2,:)-1)*sizA(1) +...
        (Im(:,3,:)-1)*sizA(1)*sizA(2);
    I0 = squeeze(I0);
    
    if existNaNmedian
        B(nchunk(i):nchunk(i+1)) = nanmedian(A(I0),2);
    else
        B(nchunk(i):nchunk(i+1)) = median(A(I0),2);
    end
end
B = cast(B,class0);
    

% % IV. DESPOT2-FM C-Function Wrapper
%   function res = despot2fm_model(x)
%     
%     % Setup FV Call
%     if fix_flag == 1
%       
%     elseif fix_flag == 0
%       fv = [x(1) 1./vox_r1 1./x(:,2) x(:,3)];
%     elseif fix_flag == 2
%       fv = [x(1) 1./vox_r1 1./x(:,2) vox_omega];
%     end
%     
%     % Call C-Code
%     res =       cpDESPOT2_residuals_SAH(fv, 0,   vox_data_0',   vox_alpha', tr, 1)';
%     res = res + cpDESPOT2_residuals_SAH(fv, 180, vox_data_180', vox_alpha', tr, 1)';
% 
%   end
% 
% % IV. DESPOT2-FM DEBUG M-Function
%   function res = despot2fm_dbg(x, dbgmode) %#ok<INUSD>
%     
%     % Pull out model parameters
%     pd_mod    = x(1);
%     r2_mod    = x(2);
%     omega_mod = x(3);
%     
%     % Flip angle terms
%     sina  = sin(vox_alpha);
%     cosa  = cos(vox_alpha);
% 
%     E1   = exp(-tr.*vox_r1);
%     E2   = exp(-tr.*r2_mod);
%     
%     % SSFP-180 Signal
%     beta = omega_mod*2*pi*tr + pi;
%     sinb = sin(beta);
%     cosb = cos(beta);
%     
%     denom = (1-E1.*cosa) .* (1-E2.*cosb) - E2.*(E1-cosa).*(E2-cosb);
%     
%     Mx = pd_mod .* ((1-E1).*E2.*sina.*sinb)      ./ denom;
%     My = pd_mod .* ((1-E1).*E2.*sina.*(cosb-E2)) ./ denom;
%     
%     ssfp_180 = sqrt(Mx.^2 + My.^2) * sqrt(exp(tr.*r2_mod));
%     
%     % SSFP-0 Signal
%     beta = omega_mod*2*pi*tr;
%     sinb = sin(beta);
%     cosb = cos(beta);
%     
%     denom = (1-E1.*cosa) .* (1-E2.*cosb) - E2.*(E1-cosa).*(E2-cosb);
%     
%     Mx = pd_mod .* ((1-E1).*E2.*sina.*sinb)      ./ denom;
%     My = pd_mod .* ((1-E1).*E2.*sina.*(cosb-E2)) ./ denom;
%     
%     ssfp_0   = sqrt(Mx.^2 + My.^2) * sqrt(exp(tr.*r2_mod));
%     
%     % Output Residual
%     if ~exist('dbgmode', 'var')
%       res = norm([ssfp_180 ssfp_0] - [vox_data_180 vox_data_0]); % SOS Res
%     else
%       res = [ssfp_180 ssfp_0];
%     end
%     
%   end
% 
% end

