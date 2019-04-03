%
%   test_simulated_world_002.m
%
%   Martin Krucinski
%
%   Evaluate suction cup model on simulated world objects & bins
%
%   2019-04-01      v001    inital version
%
%   2019-04-03      v002    using camera matrices from Corke's book

mm  = 1e-3;
Lx  = -50.0*mm;      % [mm]  Camera location
Lx=0    % Does this correspond to the images? They look very symmetric, both along x and y
Ly  = 0.0*mm;
Lz  = -1460.0*mm;
L   = [ Lx ; Ly ; Lz ];

%----------------------------------------------------------------------
%   Define grasp quality evaluation parameters
%
Std_Bad      = 0.150; %***0.030;           % standard deviation for pick points to be flagged as BAD

home_dir    = 'C:\Users\USMAKRU\OneDrive - ABB\2018\Logistics\Suction Cup Modeling\Suction_1.0_ML';
cd(home_dir)
%home_dir    = pwd;

cd('Simulated_World\Bins_001')

img_raw_color   = imread('s00000050.png');%'b00000002.png');
img_raw_depth   = imread('s00000050Depth.png');%'b00000002d.png');


%----------------------------------------------------------------------
%   Convert depth image

depth_scale     = 22e-6%* 2;     % Note! Simulated object size is 2x real object size
%***inputDepth      = -double(img_raw_depth)*depth_scale + Lz;
%inputDepth      = +double(img_raw_depth)*depth_scale - Lz;
%inputDepth      = +double(img_raw_depth)*depth_scale + Lz;

% disp('Using custom depth map conversion: 0.8 - img_raw*depth_scale')
% inputDepth      = 0.8 - double(img_raw_depth)*depth_scale;

% DOESN'T WORK inputDepth      = -double(img_raw_depth)*depth_scale - Lz;

inputDepth      = +double(img_raw_depth)*depth_scale - Lz/2;

%----------------------------------------------------------------------
%  Plot the loaded image data


N_metrics           = 2;    % number of metrics to record for each pick point

Suction_Cup_Check = true;

%----------------------------------------------------------------------
f1=figure;
image(img_raw_color);
title('Raw Color Image')
axis equal

%----------------------------------------------------------------------
f2=figure;
max_depth = 2^16;
image(double(img_raw_depth)/max_depth * 64);
colormap('gray')
title('Raw Depth Image')
axis equal

%----------------------------------------------------------------------

cd(home_dir)

%cameraIntrinsics    = load('simulated_world-intrinsics.txt','ascii');
%
%  6.17043335e+02	  0.00000000e+00	  2.96958984e+02	
%  0.00000000e+00	  6.16720093e+02	  2.39890137e+02	
%  0.00000000e+00	  0.00000000e+00	  1.00000000e+00
 
%   Construct camera intrinsics matrix directly

% [ fx 0 0 ; s fy 0 ; cx cy 1 ]

s   =   0;
f   = 20*mm %80 %20;    % [mm] camera focal length

M1 = [f 0 0 ; 0 f 0 ; 0 0 1];
M2 = [1 0 0 0 ; 0 1 0 0 ; 0 0 1 0];
M=M1*M2;

%   Sensor size

Sx  = 23.66*mm;
Sy  = 13.3*mm;

Iheight     = 720;
Iwidth      = 1280;
cx  =   Iwidth / 2;
cy  =   Iheight / 2;

% *** NOTE! The PrimeSense settings are probably IRRELEVANT for Frol's
% images
% px  = 1.4e-3;   % [mm]  PrimeSense D415 pixel size 1.4 um
% py  = 1.4e-3;   % [mm]  PrimeSense D415 pixel size 1.4 um
px  = Sx / Iwidth;   % [mm]  Frol's camera sensor
py  = Sy / Iheight;   % [mm]  Frol's camera sensor

fx  = f/px;
fy  = f/py;

%K   = [ fx 0 0 ; s fy 0 ; cx cy 1 ];
K   = [ fx s cx ; 0 fy cy ; 0 0 1 ];

CIM = K * M2;


% cameraPose          = load('test-camera-pose.txt','ascii');
% 
% 0.9999    0.0115    0.0050    0.1189
% 0.0062   -0.8006    0.5991   -0.1604
% 0.0109   -0.5990   -0.8006    0.7012
% 0         0         0         1.0000

%cameraPose          = load('simulated_world-camera-pose.txt','ascii');

%   Quaterion format of camera extrinsic matrix

s2  = 1 / sqrt(2);
Q1  = quaternion(s2,0,s2,0);
%R1  = rotmat(Q1,'frame');

%R1  = [-1 0 0 ; 0 1 0 ; 0 0 +1]';
%R1  = eye(3);

%R1      = [ 0 1 0 ; 1 0 0 ; 0 0 -1]';   % Frol's camera orientation
R1      = [ 0 1 0 ; 1 0 0 ; 0 0 -1];   % Frol's camera orientation

%cameraPose  = [ R1 L ; 0 0 0 1];
CEM  = [ R1 -L ; 0 0 0 1];

%----------------------------------------------------------------------

[worldX,worldY, worldZ] = project_D_to_XYZ(inputDepth,CIM,CEM);

%----------------------------------------------------------------------
f5=figure;
cld_world   = plot3(worldX, worldY, worldZ, 'k.');

colormap('gray')

%xlim([-1 1])
%ylim([-1 1])
%zlim([-0.5 0.5])
axis equal
xlabel('x')
ylabel('y')
title('World 3D Point Cloud')

