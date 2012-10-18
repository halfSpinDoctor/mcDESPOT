% FUNCTION [] = animateMag(M1, t, tstart)
%
% Function to animate evolution of the magnitization vector
%
% Samuel A. Hurley
% v1.0 23-Nov-2010
% University of Wisconsin
%
% Changelog:
%     v1.0 - Initial version

function [] = animateMagTwoPool(M1, t, tstart)

if ~exist('tstart', 'var')
  tstart = 1;
end

% Initiate the plot
axis;

xlabel('Mx')
ylabel('My')
zlabel('Mz')

xlim([-1 1]);
ylim([-1 1]);
zlim([-1 1]);

grid on;
axis square;

% Setup the two vectors, for faster drawing
ii = tstart;
mvec = line('XData', [0 M1(ii,1)], 'YData', [0 M1(ii,3)], 'ZData', [0 M1(ii, 5)], 'Color', 'r');
fvec = line('XData', [0 M1(ii,2)], 'YData', [0 M1(ii,4)], 'ZData', [0 M1(ii, 6)], 'Color', 'b');

legend('Myelin Water', 'Free Water');

% Animate!
for ii = tstart+1:size(M1, 1)
  set(mvec, 'XData', [0 M1(ii,1)], 'YData', [0 M1(ii,3)], 'ZData', [0 M1(ii, 5)]);
  set(fvec, 'XData', [0 M1(ii,2)], 'YData', [0 M1(ii,4)], 'ZData', [0 M1(ii, 6)]);
  pause(t);
end


return;