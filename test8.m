%
%   test8.m
%
%   Testing of Suction 1.0
%   Martin Krucinski
%
%   2019-01-23      v8

%----------------------------------------------------------------------
%   Define suction cup parameters
%
in      = 0.0254;
mm      = 1e-3;

Suction.RO      = 47/2*mm;        %   Suction cup OUTER radius
Suction.RI      = 35/2*mm;        %   Suction cup INNER radius, lip width approx. 5 mm

%----------------------------------------------------------------------
%   Define bin parameters
%

Bin.BW          = 375*mm;
Bin.BL          = 550*mm;
Bin.BH          = 220*mm;       % Bin dimensions Width x Length, Height, Origin
Bin.X           = 0.000;        % Bin origin
Bin.Y           = 0.000;
Bin.Z           = 0.000;
Bin.dx          = 0.1;          % tolerance around bin limits to use for pick point filtering
Bin.dy          = 0.1;
Bin.dz          = 0.1;

%----------------------------------------------------------------------
%   Load sensor data images for testing
%

img_raw_color   = imread('24_image_raw.color.png');
img_raw_depth   = imread('24_image_raw.depth.png');
img_proc_pick   = imread('24_image_proc.pickpoint.png');
depth_scale     = 1/10^4;
inputDepth      = double(img_raw_depth)*depth_scale;

cameraIntrinsics    = load('test-camera-intrinsics.txt','ascii');
cameraPose          = load('test-camera-pose.txt','ascii');

N_metrics           = 1;    % number of metrics to record for each pick point

Suction_Cup_Check = true;

if Suction_Cup_Check,
    disp('INCLUDING check for suction cup lip completeness...')
else
    disp('SKIPPING check for suction cup lip completeness...')
end

%----------------------------------------------------------------------
f1=figure;
image(img_raw_color);
title('Raw Color Image')

%----------------------------------------------------------------------
f2=figure;
image(img_raw_depth/256);
colormap('gray')
title('Raw Depth Image')

%----------------------------------------------------------------------
%
%   Filter depth image to move glitch points at d = 0.0 further back, to
%   e.g. 2.0 m
%
%   Filter the inputDepth image, not the img_raw_depth


d_move_back     = 20%8.0;   %   Move back distance
inputDepth_filtered     = filter_depth( inputDepth, d_move_back);

img_filtered_depth = inputDepth_filtered / depth_scale;      % create filtered depth image map

%----------------------------------------------------------------------
f3=figure;
image(img_filtered_depth/256);
colormap('gray')
title('Filtered Depth Image (glitch filtered)')

%----------------------------------------------------------------------
f4=figure;
image(img_proc_pick);
title('MIT Princeton Pick Points')

%----------------------------------------------------------------------

[worldX,worldY, worldZ] = project_D_to_XYZ(inputDepth,cameraIntrinsics,cameraPose);

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

%----------------------------------------------------------------------
%   Pickpoints for Testing with Bin sequence 24
%   0 - [0 0 2]         for testing
%   1 - [ 0.07, 0.38 0.500]   top of bottle
%   2 - [ 0.09, 0.33 0.500]   edge of bottle

PickPointT0 = [0 0 2];
PickPointT1 = [0.07 0.38 0.500];
PickPointT2 = [0.09 0.33 0.500];

PickPoints_Bin24_RAW = [...
    36.212,91.699,113.74,0,0,0,1,0.90927
    146.85,292.08,35.823,0,0,0,1,0.83563
    260.16,321.12,55.76,0,0,0,1,0.82962
    294.51,275.79,58.278,0,0,0,1,0.80612
    312.68,493.64,107.75,0,0,0,1,0.77827
    233.77,87.843,75.873,0,0,0,1,0.73422
    236.07,87.393,74.302,0,0,0,1,0.73406
    323.06,378.72,88.017,0,0,0,1,0.72525
    71.451,486.57,34.76,0,0,0,1,0.69868
    93.827,347.08,70.401,0,0,0,1,0.65517
    106.52,72.222,62.138,0,0,0,1,0.63758
    138.92,477.18,78.111,0,0,0,1,0.55286
    105.98,200.38,63.111,0,0,0,1,0.46411
    ];

PickPoints   = PickPoints_Bin24_RAW(:, 1:3)*mm;

%----------------------------------------------------------------------
%   Generate pick points for heat map

NskipX          = 10;   % pixel spacing between pick points to evaluate
NskipY          = 10;

[PickPoints,pick_set_X,pick_set_Y]    = Generate_Pick_Points(NskipX, NskipY, worldX, worldY, worldZ);


%----------------------------------------------------------------------
%   Filter out pick points that are outside bin


%***[PickPoints_Bin,pick_set_X_Bin,pick_set_Y_Bin]    = Filter_Bin_Pick_Points(PickPoints, Bin);

[PickPoints_Bin,pick_rows]    = Filter_Bin_Pick_Points(PickPoints, Bin);


%----------------------------------------------------------------------
f6=figure;

PickPoints_Bin_X        = PickPoints_Bin(:,1);
PickPoints_Bin_Y        = PickPoints_Bin(:,2);
PickPoints_Bin_Z        = PickPoints_Bin(:,3);

cld_Pick_Bin            = plot3(PickPoints_Bin_X, PickPoints_Bin_Y, PickPoints_Bin_Z, 'b.');
axis equal
colormap('gray')

X_min               = Bin.X - Bin.dx;
X_max               = Bin.X + Bin.BW + Bin.dx;

Y_min               = Bin.Y - Bin.dy;
Y_max               = Bin.Y + Bin.BL + Bin.dy;

Z_min               = Bin.Z - Bin.dz;
Z_max               = Bin.Z + Bin.BH + Bin.dz;



X_min                 = Bin.X - Bin.dx;
X_max                 = Bin.X + Bin.BW + Bin.dx;

Y_min                 = Bin.Y - Bin.dy;
Y_max                 = Bin.Y + Bin.BL + Bin.dy;

Z_min                 = Bin.Z - Bin.dz;
Z_max                 = Bin.Z + Bin.BH + Bin.dz;


X_range = [ X_min X_max ];
Y_range = [ Y_min Y_max ];
Z_range = [ Z_min Z_max ];

xlim(X_range)
ylim(Y_range)
zlim(Z_range)

xlabel('x')
ylabel('y')
title('Pick Points 3D Point Cloud')


%----------------------------------------------------------------------
%
%   Process a list of pick points
%

N_pick          = length(PickPoints_Bin);
PickMetrics     = NaN(N_pick, N_metrics);   %   Matrix to store all evaluation metrics

for pick_point_no = 1:N_pick,
    if mod(pick_point_no, 100) == 0,
        % home
        disp([ 'Processing pick point ' num2str(pick_point_no) '...' ]);
    end
    
    
    PickPoint   = PickPoints_Bin(pick_point_no,:);
    
    PickPose = [
        1 0 0 -PickPoint(1) ;
        0 1 0 -PickPoint(2) ;
        0 0 1 -PickPoint(3) ;
        0 0 0 1];   % pick point location and approach vector in T matrix form
    % currently approach vector is hard coded to be from straight above
    % along the z-direction
    
    %   Transform point cloud as viewed along the approach vector direction
    %   at the pick point location, i.e. a relative pick point point cloud
    
    
    [ Pick_PC_Rel_X, Pick_PC_Rel_Y, Pick_PC_Rel_Z ] = Transform_PC(WorldX, WorldY, WorldZ, PickPose);
    
    %   Select the point cloud that corresponds to the projection of the
    %   suction cup lip
    %   i.e. the set of 3D point cloud points whose pickX & pickY lie within
    %   the lip area
    
    [ lipX, lipY, lipZ ]            = Extract_Lip( Pick_PC_Rel_X, Pick_PC_Rel_Y, Pick_PC_Rel_Z, Suction);
    
    Std1=std(lipZ);        % Metric 1 - total Z-height standard deviation
    
    %   Calculate points for the suction cup lip circles with radii RI & RO
    
    [CIx, CIy]  = circle(Suction.RI);
    [COx, COy]  = circle(Suction.RO);
    
    %   Centerline circle is CM
    
    Suction.RM = mean([ Suction.RI Suction.RO]);    % mean suction cup lib radius
    
    [CMx, CMy]  = circle(Suction.RM);
    Lip_Complete_Flag = true;
    %*** replace wtih Lip_Complete_Flag = Lip_Complete_Check( lipX, lipY,
    %lipZ, CMx, CMy, Rcheck)
    
    
