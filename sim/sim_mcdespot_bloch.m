% FUNCTION [s_spgr s_ssfp_0 s_ssfp_180] = sim_mcdespot_bloch(fv, alpha, tr_spgr, tr_ssfp, ss_tol)
%
% Compute theoretical signal levels directly from Bloch equation numerical
% solution.  Evaluates two-pool Bloch equations with exchange (but no MT
% effect) for prototypical SPGR (Mxy re-set to 0 at end of each TR, TE=0ms) and
% SSFP (TE = TR/2, instantaneous RF pulse) sequences.
%
% Loops over multiple TRs until the change in signal is <.001%, in order to
% reach steady-state.
%
% T1m = T1 fast/"myelin"   T1f = T1 slow/"free"
%
% Inputs:
%    fv               - Parameter Vector [t1_m t1_f t2_m t2_f MWF tau omega pd_spgr pd_ssfp]
%                       t1, t2, tau in seconds, 0 < MWF < 1, omega in Hz, pd in a.u.
%    alpha            - flip angles
%    tr_spgr, tr_ssfp - repitition times
%    ss_tol           - tolerance in signal change between TR(n) and TR(n-1) to determine
%                       steady-state condition
%
% Outputs:
%    s_spgr
%    s_ssfp_0
%    s_ssfp_180
%
% Samuel A. Hurley
% University of Wisconsin
% v4.0 23-Nov-2010
%
% Changelog:
%       v4.0 Initial version, based off of sim_mcdespot_m.m

function [s_spgr s_ssfp_0 s_ssfp_180 M_spgr M_ssfp_0 M_ssfp_180] = sim_mcdespot_bloch(fv, alpha, tr_spgr, tr_ssfp, ss_tol)

% Set ODE tolerance
opts = odeset('RelTol', ss_tol);

% Convert flip angles to radians
alpha = alpha .* pi/180;

MWF = fv(5);

% Bloch Equation Function Handle
bloch = @(t,x) blochMcConnell(x, fv, 0);
M_spgr     = [];
M_ssfp_0   = [];
M_ssfp_180 = [];

tic;

%% == SPGR Sequence ==

fprintf('Working on SPGR:');

% Initalize output vector for SPGR
s_spgr = zeros(size(alpha));

% Loop over flip angles
for ii = 1:length(alpha)
  progressbar(ii/length(alpha));
  
  % Compute Rx rotation matrix
  R = [1    0    0                 0               0               0
       0    1    0                 0               0               0
       0    0    cos(alpha(ii))    0              -sin(alpha(ii))  0
       0    0    0                 cos(alpha(ii))  0              -sin(alpha(ii))
       0    0    sin(alpha(ii))    0               cos(alpha(ii))  0
       0    0    0                 sin(alpha(ii))  0               cos(alpha(ii))];
  
  % Initial Magnitization Vector
  M = [0 0 0 0 MWF (1-MWF)]' * fv(8);
  
  % Initial observed signal
  S = 1;
  S_Prev = 0;
  
  % DEBUG: Print out steady-state info from previous iter
%   if exist('n', 'var')
%     fprintf([' FA ' num2str(alpha(ii-1)*180/pi, '%02.0f') ' steady-state in ' num2str(n, '%04.0f') ' TRs']);
%   end
  
  n = 0;
  while abs(S - S_Prev)/S > ss_tol
    
    S_Prev = S;
    
    % Instantaeous Flip
    M = R*M;
    
    % Signal is abs(Mxy) at TE = 0 (ideal SPGR case)
    S = abs(M(1) + 1i*M(3)) + abs(M(2) + 1i*M(4));
    
    % Free precession over 1 TR (yes, SPGR is a type of SSFP!)
    [T M2] = ode45(bloch, [0 tr_spgr], M, opts);

    if nargout > 3
      % DEBUG Save magnitization history
      M_spgr = [M_spgr; M2]; %#ok<AGROW>
    end
    
    M = M2(end, :)';
    
    % Full spoiling -- set all transverse mag to zero
    M(1:4) = [0 0 0 0]';
    
    n = n+1;
  
  end
  
  s_spgr(ii) = S;
  
end

toc;

%% == SSFP-0 Sequence ==

fprintf('Working on SSFP-0:');

% Initalize output vector for SSFP
s_ssfp_0 = zeros(size(alpha));
clear n;

% Loop over flip angles
for ii = 1:length(alpha)
  progressbar(ii/length(alpha));
  
  % Reset Initial Magnitization Vector to M0
  M = [0 0 0 0 MWF (1-MWF)]' * fv(9);
  
  % Initial observed signal
  S = 1;
  S_Prev = 0;
  
  % DEBUG: Print out steady-state info from previous iter
%   if exist('n', 'var')
%     fprintf([' FA ' num2str(alpha(ii-1)*180/pi, '%02.0f') ' steady-state in ' num2str(n, '%04.0f') ' TRs']);
%   end
  
  n = 0;
  while abs(S - S_Prev)/S > ss_tol

    % Instantaneous Flip -- alpha/2 for first pulse
    if n == 0
      a = alpha(ii) / 2;
    else
      a = alpha(ii);
    end
    
    % Compute Rx rotation matrix
    R = [1    0    0         0        0       0
         0    1    0         0        0       0
         0    0    cos(a)    0       -sin(a)  0
         0    0    0         cos(a)   0      -sin(a)
         0    0    sin(a)    0        cos(a)  0
         0    0    0         sin(a)   0       cos(a)];
       
    M = R*M;
    
    % Free Precession for TR/2
    [T M2] = ode45(bloch, [0 (tr_ssfp/2)], M, opts);
    
    if nargout > 3
      % DEBUG Save magnitization history
      M_ssfp_0 = [M_ssfp_0; M2]; %#ok<AGROW>
    end
    
    M = M2(end, :)';
    
    % Free precession for 2nd TR/2
    [T M2] = ode45(bloch, [(tr_ssfp/2) tr_ssfp], M, opts);
    
    if nargout > 3
      % DEBUG Save magnitization history
      M_ssfp_0 = [M_ssfp_0; M2]; %#ok<AGROW>
    end
    
    M = M2(end, :)';
    
    % !! READOUT AT END OF TR !!
    % Signal is abs(Mxy) at TE = TR
    S_Prev = S;
    S = abs(M(1) + M(2) + 1i*M(3) + 1i*M(4));
    
    n = n+1;
  
  end
  
  s_ssfp_0(ii) = S;
  
end

toc;


%% == SSFP-180 Sequence ==

fprintf('Working on SSFP-180:');

% Initalize output vector for SSFP
s_ssfp_180 = zeros(size(alpha));
clear n;

% Loop over flip angles
for ii = 1:length(alpha)
  progressbar(ii/length(alpha));
  
  % Reset Initial Magnitization Vector to M0
  M = [0 0 0 0 MWF (1-MWF)]' * fv(9);
  
  % Initial observed signal
  S = 1;
  S_Prev = 0;
  
  % DEBUG: Print out steady-state info from previous iter
%   if exist('n', 'var')
%     fprintf([' FA ' num2str(alpha(ii-1)*180/pi, '%02.0f') ' steady-state in ' num2str(n, '%04.0f') ' TRs']);
%   end
  
  n = 0;
  while abs(S - S_Prev)/S > ss_tol

    % Instantaneous Flip -- alpha/2 for first pulse
    % Chop RF for SSFP-180
    if n == 0
      a = alpha(ii) / 2;
    elseif round(n/2) ~= n/2 % Odd
      a = alpha(ii);
    else
      a = -alpha(ii);        % Even
    end
    
    % Compute Rx rotation matrix
    R = [1    0    0         0        0       0
         0    1    0         0        0       0
         0    0    cos(a)    0       -sin(a)  0
         0    0    0         cos(a)   0      -sin(a)
         0    0    sin(a)    0        cos(a)  0
         0    0    0         sin(a)   0       cos(a)];
       
    M = R*M;
    
    % Free Precession for TR/2
    [T M2] = ode45(bloch, [0 (tr_ssfp/2)], M, opts);
    
    if nargout > 3
      % DEBUG Save magnitization history
      M_ssfp_180 = [M_ssfp_180; M2]; %#ok<AGROW>
    end
    
    M = M2(end, :)';
    
    % Free precession for 2nd TR/2
    [T M2] = ode45(bloch, [(tr_ssfp/2) tr_ssfp], M, opts);
    
    if nargout > 3
      % DEBUG Save magnitization history
      M_ssfp_180 = [M_ssfp_180; M2]; %#ok<AGROW>
    end
    
    M = M2(end, :)';
    
    % !! Readout at the end of the TR !!
    % Signal is abs(Mxy) at TE = TR
    S_Prev = S;
    S = abs(M(1) + M(2) + 1i*M(3) + 1i*M(4));
    
    n = n+1;
  
  end
  
  s_ssfp_180(ii) = S;
  
end

toc;
