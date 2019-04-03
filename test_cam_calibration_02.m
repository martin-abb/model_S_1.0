%   test_cam_calibration_02.m
%
%   run after CEM (camera extrinsics) and CIM (camera intrinsics) have been
%   defined
%
%   2019-04-03      

%W1  = [0 0 0]'     % origin

%W1  = [+50*mm 0 0]'     % origin
%W1rel   = [0.1 0.1 0]'
%W1rel     = [ +0.3 +0.2 +0.0]'
%W1rel   = [0.6 0 -0.3]'

W1rel   = [0.3142 -0.5595 +0.0205]'         % bottom of bin corner
%W1rel   = [0.2542 -0.5395 +0.4580]'         % top of bin corner

%W1  = [-50*mm 0 0]' + W1rel     % absolute world coordinate
W1  = [0 0 0]' + W1rel     % absolute world coordinate

W1x = [ W1 ; 1];    % extended

%--------------------------------------------------------------------
%   Per Corke's definition

% W2      = [ 0.3 0.4 3.0 ]';
% W2x     = [ W2 ; 1];
% 
% ptemp1  = M * W2x;
% ptilde1 = ptemp1(1:2) / ptemp1(3)


% L2      = [ -0.5 0 0 ]';
% M2t    = [ eye(3) -L2 ];
% Mb      = M1 * M2t ;
% ptemp2  = Mb * W2x;
% ptilde2 = ptemp2(1:2) / ptemp2(3)
% 


%--------------------------------------------------------------------
%   Previous definitions

C       = CIM * CEM ; 
ptilde  = C*W1x;

Z       = ptilde(3);
p       = round( ptilde / Z );

pixx      = p(1);
pixy      = p(2);


disp('Pixel coordinates')
disp([ 'pixx  = ' num2str(pixx) ', pixy = ' num2str(pixy) ', Z = ' num2str(Z) ])
pixy_I = Iheight - pixy;
disp([ 'Image y-coord = Iheight - pixy = ' num2str(pixy_I) ])


