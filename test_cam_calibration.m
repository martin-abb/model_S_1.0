%   test_cam_calibration.m
%
%   run after CEM (camera extrinsics) and CIM (camera intrinsics) have been
%   defined

%W1  = [0 0 0]'     % origin

%W1  = [+50*mm 0 0]'     % origin
W1  = [+50*mm 0 0]' + [ -0.3 0.2 0.3]'     % origin


W1x = [ W1 ; 1]    % extended
C1x = CEM*W1x
C1  = C1x(1:3)

% T1  = CIM'*C1
% Z   = T1(3)
% 
% P1  = [ T1(1)/Z T1(2)/Z Z]
%T1  = CIM'*C1

Z   = C1(3)
%T1  = [ C1(1)/Z C1(2)/Z 1 ]'
T1  = C1 / Z

P1  = CIM' * T1
P1(3)       = Z

disp('Pixel coordinates')
disp([ 'Px  = ' num2str(P1(1)) ', Py = ' num2str(P1(2)) ', Pz = ' num2str(P1(3)) ])

