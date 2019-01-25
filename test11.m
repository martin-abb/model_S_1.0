%
%   test11.m
%
%   Testing of Suction 1.0
%   Martin Krucinski
%
%   2019-01-23      v8
%   2019-01-24      v9  Added subtraction of point clouds for pick point
%                       lips to increase their visibility
%   2019-01-25      v10 Speed ups for dense heat map generation
%   2019-01-25      v11 Include frequency / slope effect of the tangential
%                       direction of the suction cup lip

%   Define running mode

script_flag = 0;        % 1 - running as multiple pick points evaluation script, DO NOT plot individual point cloud plots
% 0 - DO plot individual point cloud plots

plot_PC_overlay = 1;    % 1 - DO plot point cloud with suction cup lip points overlaid
% 0 - DO NOT plot the enhanced point cloud

NskipX          = 200%10;   % pixel spacing between pick points to evaluate
NskipY          = 200%10;

%----------------------------------------------------------------------
%   Define suction cup parameters
%
in      = 0.0254;
mm      = 1e-3;

Suction.RO      = 20*mm %*****47/2*mm;        %   Suction cup OUTER radius
Suction.RI      = 15*mm %***35/2*mm;        %   Suction cup INNER radius, lip width approx. 5 mm

%----------------------------------------------------------------------
%   Define bin parameters
%

Bin.BW          = 375*mm;
Bin.BL          = 550*mm;
Bin.BH          = 220*mm;       % Bin dimensions Width x Length, Height, Origin
Bin.X           = 0.000;        % Bin origin
Bin.Y           = 0.000;
Bin.Z           = 0.000;

% Bin.dx          = 0.1;          % tolerance around bin limits to use for pick point filtering
% Bin.dy          = 0.1;
% Bin.dz          = 0.1;

Bin.dx          = Suction.RO;          % tolerance around bin limits to use for pick point filtering
Bin.dy          = Suction.RO;
Bin.dz          = 0.01;

%   The ranges below were used for point cloud axis limits, extend BEYOND
%   bin
% Bin.X_min                 = Bin.X - Bin.dx;
% Bin.X_max                 = Bin.X + Bin.BW + Bin.dx;
%
% Bin.Y_min                 = Bin.Y - Bin.dy;
% Bin.Y_max                 = Bin.Y + Bin.BL + Bin.dy;
%
% Bin.Z_min                 = Bin.Z - Bin.dz;
% Bin.Z_max                 = Bin.Z + Bin.BH + Bin.dz;

%   Use these ranges for checking if pick point inside bin, REDUCE the bin
%   size instead of extending

Bin.X_min                 = Bin.X + Bin.dx;
Bin.X_max                 = Bin.X + Bin.BW - Bin.dx;

Bin.Y_min                 = Bin.Y + Bin.dy;
Bin.Y_max                 = Bin.Y + Bin.BL - Bin.dy;

Bin.Z_min                 = Bin.Z + Bin.dz;
Bin.Z_max                 = Bin.Z + Bin.BH - Bin.dz;

%----------------------------------------------------------------------
%   Define grasp quality evaluation parameters
%
Std_Bad      = 0.150; %***0.030;           % standard deviation for pick points to be flagged as BAD

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

N_metrics           = 2;    % number of metrics to record for each pick point

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

[PickPoints,pick_set_X,pick_set_Y]    = Generate_Pick_Points(NskipX, NskipY, worldX, worldY, worldZ);


%----------------------------------------------------------------------
%   Filter out pick points that are outside bin


%***[PickPoints_Bin,pick_set_X_Bin,pick_set_Y_Bin]    = Filter_Bin_Pick_Points(PickPoints, Bin);

if 0, %*** CAN'T FILTER IN ORDER TO GET A RECTANGULAR HEAT MAP BACK, SEE COMMENTS BELOW
    [PickPoints_Bin,pick_rows]    = Filter_Bin_Pick_Points(PickPoints, Bin);
else
    PickPoints_Bin          = PickPoints;
end


%----------------------------------------------------------------------
f6=figure;

PickPoints_Bin_X        = PickPoints_Bin(:,1);
PickPoints_Bin_Y        = PickPoints_Bin(:,2);
PickPoints_Bin_Z        = PickPoints_Bin(:,3);

cld_Pick_Bin            = plot3(PickPoints_Bin_X, PickPoints_Bin_Y, PickPoints_Bin_Z, 'b.');
axis equal
colormap('gray')

plot_X_range = [ Bin.X - 0.1 Bin.X + Bin.BW + 0.1 ];
plot_Y_range = [ Bin.Y - 0.1 Bin.Y + Bin.BL + 0.1 ];
plot_Z_range = [ Bin.Z - 0.1 Bin.Z + Bin.BH + 0.1 ];

xlim(plot_X_range)
ylim(plot_Y_range)
zlim(plot_Z_range)

xlabel('x')
ylabel('y')
title('Pick Points 3D Point Cloud')


%----------------------------------------------------------------------
%
%   Process a list of pick points
%

all_linear_Ind_Union  = [];
all_linear_Ind  = {};
all_lipX        = {};
all_lipY        = {};
all_lipZ        = {};

N_pick          = length(PickPoints_Bin);
PickMetrics     = NaN(N_pick, N_metrics);   %   Matrix to store all evaluation metrics

all_Lip_Complete_Flag = zeros(N_pick,1);

disp([ 'Total number of pick points: ' num2str(N_pick) ]);

for pick_point_no = 1:N_pick,
    if mod(pick_point_no, 100) == 0,
        % home
        disp([ 'Processing pick point ' num2str(pick_point_no) '...' ]);
    end
    
    
    PickPoint   = PickPoints_Bin(pick_point_no,:);
    PX = PickPoint(1);
    PY = PickPoint(2);
    PZ = PickPoint(3);
    
    PickPose = [
        1 0 0 -PX ;
        0 1 0 -PY ;
        0 0 1 -PZ;
        0 0 0 1];   % pick point location and approach vector in T matrix form
    % currently approach vector is hard coded to be from straight above
    % along the z-direction
    
    %   Transform point cloud as viewed along the approach vector direction
    %   at the pick point location, i.e. a relative pick point point cloud
    
    
    [ Pick_PC_Rel_X, Pick_PC_Rel_Y, Pick_PC_Rel_Z ] = Transform_PC(worldX, worldY, worldZ, PickPose);
    
    %   Select the point cloud that corresponds to the projection of the
    %   suction cup lip
    %   i.e. the set of 3D point cloud points whose pickX & pickY lie within
    %   the lip area
    
    [ lipX, lipY, lipZ, linearInd ]            = Extract_Lip2(Pick_PC_Rel_X, Pick_PC_Rel_Y, Pick_PC_Rel_Z, Suction);
    
    all_linear_Ind_Union            = union(all_linear_Ind_Union, linearInd);
    all_linear_Ind{pick_point_no}   = linearInd;
    
    %   convert stored all_lip_X etc. to ABSOLUTE coordinates, not RELATIVE
    all_lip_X{pick_point_no}        = lipX + PX;
    all_lip_Y{pick_point_no}        = lipY + PY;
    all_lip_Z{pick_point_no}        = lipZ + PZ;
    
    
    % **** DONE ****
    %*** expand Extract_Lip function to also return the indices of
    % extracted points for subsequent processing and set removal of these
    % points from the World 3D Point Cloud for plotting purposes
    
    Std1=std(lipZ);        % Metric 1 - total Z-height standard deviation
    
    %   Calculate points for the suction cup lip circles with radii RI & RO
    
    [CIx, CIy]  = circle(Suction.RI);
    [COx, COy]  = circle(Suction.RO);
    
    %   Centerline circle is CM
    
    Suction.RM      = mean([ Suction.RI Suction.RO]);    % mean suction cup lib radius
    Suction.Rcheck  = (Suction.RO - Suction.RI)/2;      % radius distance to check around for lip points, half of the suction lip width
    
    %***    Lip_Complete_Flag = true;
    Lip_Complete_Flag = Lip_Complete_Check( lipX, lipY, lipZ, Suction);
    
    all_Lip_Complete_Flag(pick_point_no) = Lip_Complete_Flag;
    
    if ~Lip_Complete_Flag,
        Std1 = Std_Bad;
    end
    
    %   Try direct implementation instead of function to check speed
    %   166 s with Nskip = 10 (on battery power) with FUNCTION
    %   With the direct implementation, it took 165 s, NO DIFFERENCE,
    %   revert back to function implementation
    
    if 1,
        [PickPoint_Bin,rows]  = Filter_Bin_Pick_Points(PickPoint, Bin);
        
        %----------------------------------------------------------------------
        %   Check if pick point is within bin limits
        
        if isempty(PickPoint_Bin),
            Std1 = Std_Bad;
        end
    else
        
        %   PickPoint outside bin
        
        if ~( (PX >= Bin.X_min) & (PX <= Bin.X_max) & ...
                (PY >= Bin.Y_min) & (PY <= Bin.Y_max) & ...
                (PZ >= Bin.Z_min) & (PZ <= Bin.Z_max) ),
            Std1 = Std_Bad;
        end
        
    end
    
    
    
    %----------------------------------------------------------------------
    %   Sort lipZ values along lip angle
    
    [lipZSorted, lipAngles, lipZ_sort_indices ] = Sort_Lip_Z(lipX, lipY, lipZ, Suction);

    
    %----------------------------------------------------------------------
    %  Include frequency / slope effect of the tangential
    %  direction of the suction cup lip
    
    
    
    
    %----------------------------------------------------------------------
    
    if ~script_flag,
        f7      = figure;
        plot(lipZSorted)
        title(['Suction cup lip Z-heights for Pick Point ' num2str(pick_point_no) ])
        xlabel('lip angle sequence #')
    end
    
    %----------------------------------------------------------------------
    
    
    if ~script_flag,
        f7b      = figure;
        plot(lipAngles, lipZSorted)
        title(['Suction cup lip Z-heights for Pick Point ' num2str(pick_point_no) ])
        xlabel('lip angle [rad]')
    end
    
    %----------------------------------------------------------------------   
    if ~script_flag,
        f8      = figure;
        
        subplot(221)
        plot3(lipX + PickPoint(1), lipY + PickPoint(2), lipZ + PickPoint(3), 'r.');
        t81=text(PickPoint(1),PickPoint(2),PickPoint(3),num2str(pick_point_no));
        xlabel('x')
        ylabel('y')
        view(3)
        title([ 'Pick Point number ' num2str(pick_point_no) ])
        
        subplot(222)
        plot3(lipX + PickPoint(1), lipY + PickPoint(2), lipZ + PickPoint(3), 'r.');
        t82=text(PickPoint(1),PickPoint(2),PickPoint(3),num2str(pick_point_no));
        xlabel('x')
        ylabel('y')
        view(0, 90)
        title('From top')
        
        subplot(223)
        plot3(lipX + PickPoint(1), lipY + PickPoint(2), lipZ + PickPoint(3), 'r.');
        t83=text(PickPoint(1),PickPoint(2),PickPoint(3),num2str(pick_point_no));
        xlabel('x')
        ylabel('y')
        view(90,0)
        title('From side')
        
        subplot(224)
        t84=text(PickPoint(1),PickPoint(2),PickPoint(3),num2str(pick_point_no));
        plot3(lipX + PickPoint(1), lipY + PickPoint(2), lipZ + PickPoint(3), 'r.');
        xlabel('x')
        ylabel('y')
        view(0,0)
        title('From front')
    end
    
    if plot_PC_overlay,
        figure(f5);
        hold on
        plot3(lipX + PickPoint(1), lipY + PickPoint(2), lipZ + PickPoint(3), 'r.');
        t51=text(PickPoint(1),PickPoint(2),PickPoint(3),num2str(pick_point_no));
    end
    
    %----------------------------------------------------------------------
    %   Save metrics for this pick point
    
    PickMetrics(pick_point_no, 1)       = Std1;
end

%----------------------------------------------------------------------
%   Process metrics

%   Normalized metric: M1_scaled (0 - Bad, 1 - Good) scaled with respect to
%   max value of M1
NX          = length(pick_set_X);
NY          = length(pick_set_Y);
M1          = PickMetrics(:,1);
M1_min      = min(M1);
M1_max      = max(M1);
M1_delta    = M1_max - M1_min;

M1_scaled   = (1 - (M1 - M1_min)/M1_delta);
Image_M1_scaled    = reshape(M1_scaled, NY, NX)*63;

M1_max_abs  = 0.050;    % Absolute, max variation that suction cup can pick up

M1_abs_scaled      = min(1.0, max((1 - (M1 - M1_min)/(M1_max_abs - M1_min)),0));
Image_M1_abs_scaled    = reshape(M1_abs_scaled, NY, NX)*63;

%*** Do note that to get a rectangular metrics image back,
%*** I CAN'T filter out pick points in the image that are OUTSIDE THE BIN,
%*** ALL pickpoints in the 2D grid HAVE TO BE USED!!
%*** That means the line
%*** [PickPoints_Bin,pick_rows]    = Filter_Bin_Pick_Points(PickPoints, Bin);
%*** above CAN'T be used!

%----------------------------------------------------------------------
%   Prepare the 3D point cloud with suction cup lip point clouds REMOVED

S               = size(worldX);
NworldX         = S(1)*S(2);
all_Ind         = 1:NworldX;

ind_non_pick_points = setdiff( all_Ind, all_linear_Ind_Union );

%----------------------------------------------------------------------
%   Plot world 3D point cloud with pick point suction cup lips point clouds highlighted
%
%   Perform this by first plotting point cloud in BLACK for all NON-LIPS
%   point cloud points
%
%   Then fill in the point cloud points in RED

if plot_PC_overlay,
    
    f9=figure;
    worldX_non_pick_points     = worldX(ind_non_pick_points);
    worldY_non_pick_points     = worldY(ind_non_pick_points);
    worldZ_non_pick_points     = worldZ(ind_non_pick_points);
    
    worldX_pick_points     = worldX(all_linear_Ind_Union);
    worldY_pick_points     = worldY(all_linear_Ind_Union);
    worldZ_pick_points     = worldZ(all_linear_Ind_Union);
    
    cld_world_non_pick_points   = plot3(worldX_non_pick_points, worldY_non_pick_points, worldZ_non_pick_points, 'k.');
    hold on
    cld_world_pick_points   = plot3(worldX_pick_points, worldY_pick_points, worldZ_pick_points, 'r.');
    hold off
    
    colormap('gray')
    
    axis equal
    xlabel('x')
    ylabel('y')
    title('World 3D Point Cloud with Suction Cup Lip Point Clouds Overlaid')
    xlim(plot_X_range)
    ylim(plot_Y_range)
    zlim(plot_Z_range)
    
end


%----------------------------------------------------------------------
%   Plot world 3D point cloud with pick point suction cup lips point clouds highlighted
%
%   Perform this by first plotting point cloud in BLACK for all NON-LIPS
%   point cloud points
%
%   Then fill in the point cloud points in RED

if plot_PC_overlay,
    
    f9b=figure;
    worldX_non_pick_points     = worldX(ind_non_pick_points);
    worldY_non_pick_points     = worldY(ind_non_pick_points);
    worldZ_non_pick_points     = worldZ(ind_non_pick_points);
    
    worldX_pick_points     = worldX(all_linear_Ind_Union);
    worldY_pick_points     = worldY(all_linear_Ind_Union);
    worldZ_pick_points     = worldZ(all_linear_Ind_Union);
    
    cld_world_non_pick_points   = plot3(worldX_non_pick_points, worldY_non_pick_points, worldZ_non_pick_points, 'k.');
    hold on
    %cld_world_pick_points   = plot3(worldX_pick_points, worldY_pick_points, worldZ_pick_points, 'r.');
    for pick_point_no = 1:N_pick,
        if all_Lip_Complete_Flag(pick_point_no),
            c_lip = 'g.';
        else
            c_lip = 'r.';
        end
        
        h_lip = plot3( all_lip_X{pick_point_no}, all_lip_Y{pick_point_no}, all_lip_Z{pick_point_no}, c_lip);
       
        t91=text(PickPoints_Bin(pick_point_no,1), ...
            PickPoints_Bin(pick_point_no,2), ...
            PickPoints_Bin(pick_point_no,3), ...
            [ num2str(pick_point_no) ' : ' num2str(M1(pick_point_no)) ] );
        set(t91,'FontWeight','Bold');
        set(t91,'Color','b');
    end
    hold off
    
    colormap('gray')
    
    axis equal
    xlabel('x')
    ylabel('y')
    title('World 3D Point Cloud with Suction Cup Lip Point Clouds Overlaid')
    xlim(plot_X_range)
    ylim(plot_Y_range)
    zlim(plot_Z_range)
    
end


%----------------------------------------------------------------------
%   Plot RAW heat maps
f10 = figure;
image(Image_M1_scaled)
colormap hot
title('RAW Heat Map (normalized)')

%----------------------------------------------------------------------
%   Plot RAW heat maps
f11 = figure;
image(Image_M1_abs_scaled)
colormap hot
title('RAW Heat Map (abs scaled)')


%----------------------------------------------------------------------
%   Generate large heat map images, of same size as original camera images
%   so that they can be displayed overlaid on top of each other

Image_M1_scaled_Large      = Image_Block_Expand(Image_M1_scaled, NskipX, NskipY);
Image_M1_abs_scaled_Large  = Image_Block_Expand(Image_M1_abs_scaled, NskipX, NskipY);


%----------------------------------------------------------------------
%   Plot heat map overlaid on original RGB image with linear ramp heat map
%

f12 = figure;


h_f12_imraw     = image(img_raw_color);
hold on
h_f12_im1       = image(Image_M1_scaled_Large);

set(h_f12_imraw,'AlphaData', 0.5);
set(h_f12_im1,'AlphaData', 0.5);

Ntop    = 20;

map     = [ zeros(64-Ntop,1); [(1/Ntop):1/Ntop:1]']*[1 1 0];

map2    = [     ones(64-Ntop,1) * [0 0 1] ;
    [(1/Ntop):1/Ntop:1]'*[1 1 0] ; ...
    ];
colormap(map2)
title(['RAW Heat Map Overlay with LINEAR ramp color map with Ntop = ' num2str(Ntop) ])


%----------------------------------------------------------------------
%   Plot heat map overlaid on original RGB image with linear ramp heat map
%

f13 = figure;


h_f13_imraw     = image(img_raw_color);
hold on
h_f13_im1       = image(Image_M1_abs_scaled_Large);

set(h_f13_imraw,'AlphaData', 0.5);
set(h_f13_im1,'AlphaData', 0.5);

Ntop    = 20;

map     = [ zeros(64-Ntop,1); [(1/Ntop):1/Ntop:1]']*[1 1 0];

map2    = [     ones(64-Ntop,1) * [0 0 1] ;
    [(1/Ntop):1/Ntop:1]'*[1 1 0] ; ...
    ];
colormap(map2)
title(['RAW Heat Map (abs scaled) Overlay with LINEAR ramp color map with Ntop = ' num2str(Ntop) ])


%----------------------------------------------------------------------
%   Plot SMOOTHED heat map overlaid on original RGB image with linear ramp heat map
%

f14 = figure;

Nsmooth     = 5;
t1=NskipX^(1/Nsmooth);

imgA        = Image_M1_abs_scaled;
imgB        = imgA;
for s=1:Nsmooth,
    imgB=imresize(imgB,t1,'bilinear');
end
Image_M1_abs_scaled_smooth  = imgB;

h_f14_imraw     = image(img_raw_color);
hold on
h_f14_im1       = image(Image_M1_abs_scaled_smooth);

set(h_f14_imraw,'AlphaData', 0.5);
set(h_f14_im1,'AlphaData', 0.5);


colormap(map2)
title(['SMOOTHED Heat Map (abs scaled) Overlay with LINEAR ramp color map with Ntop = ' num2str(Ntop) ])
