% FUNCTION [s_spgr s_ssfp] = sim_mcdespot_bloch_single(fv, alpha, tr_spgr, tr_ssfp, ss_tol)
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
%    fv               - [t1 t2 omega pd] 
%                       t1, t2, tau in seconds, 0 < MWF < 1, omega in Hz, pd in a.u.
%    alpha            - flip angle
%    tr_spgr, tr_ssfp - repitition times
%    phase_cyc        - RF phase cycling
%    tau              - pulse width (set to 0 for instant pulse model)
%    ss_tol           - tolerance in signal change between TR(n) and TR(n-1) to determine
%                       steady-state condition
%
% Outputs:
%    s_spgr
%    s_ssfp
%
% Samuel A. Hurley
% University of Wisconsin
% v1.0 25-Mar-2014
%
% Changelog:
%       v1.0 Initial version, based off of sim_mcdespot_bloch_2.m

function [s_spgr s_ssfp M_spgr M_ssfp] = sim_mcdespot_bloch_single(fv, alpha, tr_spgr, tr_ssfp, tau, ss_tol)

% GAMMA = 42.577480; % MHz/T
GAMMA = 267.513e6; % rad/s/T   1H Gyromagnetic Ratio

if ~check_var('ss_tol')
  ss_tol = 1e-3;
end

% Set ODE tolerance
opts = odeset('RelTol', ss_tol);

% Convert flip angles to radians
alpha = alpha .* pi/180;

% Bloch Equation Function Handle
bloch = @(t,x)blochSingle(x, fv, 0);
M_spgr = [];
M_ssfp = [];

tic;

%% == SPGR Sequence ==

fprintf('Working on SPGR:');

% Initalize output vector for SPGR
s_spgr = zeros(size(alpha));

% Loop over flip angles
for ii = 1:length(alpha)
  % progressbar(ii/length(alpha));
  
  % Instantaneous Pulse
  % Compute Rx rotation matrix
  R = [1   0                0         
       0   cos(alpha(ii))  -sin(alpha(ii))
       0   sin(alpha(ii))   cos(alpha(ii))];
  
  % Initial Magnitization Vector
  M = [0 0 1]' * fv(4);
  
  % Initial observed signal
  S = 1;
  S_Prev = 0;
  
%  % DEBUG: Print out steady-state info from previous iter
%   if exist('n', 'var')
%     fprintf([' FA ' num2str(alpha(ii-1)*180/pi, '%02.0f') ' steady-state in ' num2str(n, '%04.0f') ' TRs']);
%   end
  
  n = 0;
  while abs(S - S_Prev)/S > ss_tol
    
    S_Prev = S;
    
    % Instantaeous Flip
    M = R*M;
    
    % Signal is abs(Mxy) at TE = 0 (ideal SPGR case)
    S = norm(M(1:2));
    
    % Free precession over 1 TR (yes, SPGR is a type of SSFP!)
    [T M2] = ode45(bloch, [0 tr_spgr], M, opts);

    if nargout > 3
      % DEBUG Save magnitization history
      M_spgr = [M_spgr; M2]; %#ok<AGROW>
    end
    
    M = M2(end, :)';
    
    % Full spoiling -- set all transverse mag to zero
    M(1:2) = [0 0]';
    
    n = n+1;
  
  end
  
  s_spgr(ii) = S;
  
end

toc;

%% == SSFP Sequence ==

fprintf('Working on SSFP-0:');

% Initalize output vector for SSFP
s_ssfp = zeros(size(alpha));
clear n;

% Loop over flip angles
for ii = 1:length(alpha)
  %progressbar(ii/length(alpha));
  
  % Reset Initial Magnitization Vector to M0
  M = [0 0 1]' * fv(4);
  
  % Initial observed signal
  S = 1;
  S_Prev = 0;
  
%  % DEBUG: Print out steady-state info from previous iter
%   if exist('n', 'var')
%     fprintf([' FA ' num2str(alpha(ii-1)*180/pi, '%02.0f') ' steady-state in ' num2str(n, '%04.0f') ' TRs']);
%   end


  a = alpha(ii);
    
  if tau == 0

    % Instantaneous Pulse
    % Compute Rx rotation matrix
    R = [1   0        0
         0   cos(a)  -sin(a)
         0   sin(a)   cos(a) ];

  else
    % Finite Pulse

    % Compute B1 amplitude of RF such that flip angle ALPHA is achieved
    % alpha(ii) = 2*pi*GAMMA*tau*b1_rf;

    b1_rf = sin(a) ./ (2*pi*GAMMA*tau);

    blochRF = @(t,x)blochSingle(x, fv, b1_rf);

  end

  n = 0;
  while abs(S - S_Prev)/S > ss_tol
    
    % Apply rotation RF pulse
    if tau == 0
      M = R*M;
    else
      [T M2] = ode45(blochRF, [0 tau], M, opts);
      M = M2(end, :)';
    end
   
    
    % Free Precession for TR/2
    [T M2] = ode45(bloch, [tau (tr_ssfp/2)], M, opts);
    
    if nargout > 3
      % DEBUG Save magnitization history
      M_ssfp = [M_ssfp; M2]; %#ok<AGROW>
    end
    
    M = M2(end, :)';
    
    % Signal is abs(Mxy) at TE = TR/2 (ideal SSFP case)
    S_Prev = S;
    S = norm(M(1:2));
    
    % Free precession for 2nd TR/2
    [T M2] = ode45(bloch, [(tr_ssfp/2) tr_ssfp], M, opts);
    
    if nargout > 3
      % DEBUG Save magnitization history
      M_ssfp = [M_ssfp; M2]; %#ok<AGROW>
    end
    
    M = M2(end, :)';
    
    n = n+1;
  
  end
  
  s_ssfp(ii) = S;
  
end


toc;
