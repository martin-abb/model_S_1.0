%
%   test_Objects.m
%
%   Martin Krucinski
%

%   2019-03-13 test_Objects.m
%   2019-03-14 Improving addObject to use linear distances
%
%   Previous version 2019-02-05 test_Wedge.m
%
%
%   Script to test suction model with auto-generated object point clouds

in      = 0.0254;
mm      = 1e-3;
dx      = 1*mm;         % world grid step size
L       = 600*mm;       % world size
W       = 400*mm;

    Std_Bad      = 0.150; %***0.030;           % standard deviation for pick points to be flagged as BAD

version = 3             % wedge version 1 - one single step
% 2 - two steps
% 3 - ramp

%R       = 0.030;        % Radius
%**** DON'T USE THIS with test_Objecst and Wedge2!!
%*** All wedges end up with the SAME SIZE in "pixels" or point cloud points
%*** use dx instead!
%***N       = 50;         % number of points along x & y directions

h       = 0.050;        % cylinder base height
A       = 0.0125%0.010;        % ramp height


init_Suction

slope   = Suction.max_slope_linear / (1/2);

%[wedgeX, wedgeY, wedgeZ] = Wedge(R, A, N, slope, version);

                all_linear_Ind_Union  = [];
                all_linear_Ind  = {};
                all_lipX        = {};
                all_lipY        = {};
                all_lipZ        = {};
                
%----------------------------------------------------------------------
%   Empty world

rangeX              = (-L/2):dx:(L/2);
rangeY              = (-W/2):dx:(W/2);

                    plot_X_range = [ rangeX(1) - 0.1 , rangeX(end) + 0.1 ];
                    plot_Y_range = [ rangeY(1) - 0.1 , rangeY(end) + 0.1 ];
                    %*** Z will be defined later, after world is known   plot_Z_range = [ Bin.Z - 0.1 Bin.Z + Bin.BH + 0.1 ];

NX                  = length(rangeX);
NY                  = length(rangeY);

[emptyX,emptyY]     = meshgrid( rangeX , rangeY );
emptyZ              = emptyX .* 0.00;   % flat plane

empty.X             = emptyX;
empty.Y             = emptyY;
empty.Z             = emptyZ;
empty.rangeX        = rangeX;
empty.rangeY        = rangeY;
emtpy.dx            = dx;

%----------------------------------------------------------------------
%   Set up all world objects

midX        = NX/2;
midY        = NY/2;


slope   = Suction.max_slope_linear / (1/2);


% dimensions in [m]
% Wedge2(        R,         h,      A,      slope_linear,   version,    dx)
wedge1  = Wedge2(0.030,     0.050,  A,       slope,         version,    dx);
wedge1_insert_pos    = [-L/2 ;-W/2 ] %[100 ; 200];

wedge2  = Wedge2(0.020,     0.040,  0.005,   slope*5,       version,    dx);
wedge2_insert_pos    = [-L/2 ; 0 ] %[100 ; 200];

wedge3  = Wedge2(0.010,     0.030,  0.010,   slope/2,       version,    dx);
wedge3_insert_pos    = [0 ; -W/2 ] %[100 ; 200];

wedge4  = Wedge2(0.015,     0.020,  0.010,   slope*2,       1,          dx);
wedge4_insert_pos    = [0.08 ; 0 ] %[100 ; 200]

wedge5  = Wedge2(0.040,     0.050,  0.020,   slope,         1,          dx);
wedge5_insert_pos    = [-0.040 ; -0.040 ] %[100 ; 200]

cylinder1   = Cylinder(0.050, 0.150, dx);
cylinder1_insert_pos    =  [0.200 0];

sphere1     = Sphere( 0.050, dx);
sphere1_insert_pos     = [0.100 -0.150];

box1        = Box(0.100, 0.040, 0.030, dx);
box1_insert_pos     = [-0.200 -0.150];

world       = addObject( empty , wedge1, wedge1_insert_pos );
world       = addObject( world , wedge2, wedge2_insert_pos );
world       = addObject( world , wedge3, wedge3_insert_pos );
world       = addObject( world , wedge4, wedge4_insert_pos );
world       = addObject( world , wedge5, wedge5_insert_pos );
world       = addObject( world , cylinder1, cylinder1_insert_pos );
world       = addObject( world , sphere1,sphere1_insert_pos );
world       = addObject( world , box1,box1_insert_pos );

worldX      = world.X;
worldY      = world.Y;
worldZ      = world.Z;
%----------------------------------------------------------------------
%   Specify pick points

all_PickPoints  = [ 0.000   0.000   0.090   ; ...
                    0.015   -0.100  0.060   ; ...
                    0.100   0.010   0.060   ; ...
                        0.0953    0.0156 0.060   ; ...
    0.2526    0.0480 0.060   ; ...
    0.1560   -0.0996 0.060   ; ...
    0.0117   -0.1899 0.060   ; ...
   -0.1513   -0.1280 0.060   ; ...
   -0.2680   -0.1675 0.060   ; ...
   -0.2806    0.0215 0.060   ; ...
    0.0138    0.0343 0.060   ; ...
                  ];
              


%----------------------------------------------------------------------
%   Plot 3D point cloud

f1=figure;
cld_world   = plot3(worldX, worldY, worldZ, 'k.');

colormap('gray')

axis equal
xlabel('x')
ylabel('y')
title('World 3D Point Cloud')
hold on

%----------------------------------------------------------------------
%   Loop over list of pick points


              
                              PickPoints_Bin          = all_PickPoints;

S_pick          = size(all_PickPoints);
N_pick          = S_pick(1);            % total number of pick points

        N_metrics           = 4;    % number of metrics to record for each pick point

                PickMetrics     = NaN(N_pick, N_metrics);   %   Matrix to store all evaluation metrics
                
                all_Lip_Complete_Flag = zeros(N_pick,1);
                
                disp([ 'Total number of pick points: ' num2str(N_pick) ]);


for pick_point_no = 1:N_pick,
    
    %PickPoint   = [0 0 0];  %   Wedge pick point is at origin
    %PickPoint   = [-0.05 0.000 0.000 ];
    
    PickPoint   = all_PickPoints(pick_point_no,:);
    PickPoint_Bin   = PickPoint;
    
    %----------------------------------------------------------------------
    %   Evaluate suction grasp score
    
    
    
    
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
                    

    %----------------------------------------------------------------------
    %   Plot 3D point cloud for suction cup lip
    
    f2=figure;
    cld_lip   = plot3(lipX, lipY, lipZ, 'r.');
    
    colormap('gray')
    
    axis equal
    xlabel('x')
    ylabel('y')
    zlim([ -0.05 0.10 ])
    title('Suction Cup Lip Point Cloud')
    
    
    
    %----------------------------------------------------------------------
    
    %***    Lip_Complete_Flag = true;
    Lip_Complete_Flag = Lip_Complete_Check( lipX, lipY, lipZ, Suction);
    
                        all_Lip_Complete_Flag(pick_point_no) = Lip_Complete_Flag;

    %Score = Suction_Model_1_0( lipX, lipY, lipZ, Suction);
    
    Score_metric = Suction_Model_1_0b( lipX, lipY, lipZ, Suction);
    
    Score       = Score_metric(1);
    Score_freq  = Score_metric(2);
    Score_amp   = Score_metric(3);
    
                        Std1=std(lipZ);        % Metric 1 - total Z-height standard deviation

    if ~Lip_Complete_Flag,
        Std1 = Std_Bad;
        Score   = 0.00;
    end
    
    disp(['Final Suction Cup Grasp Score = ' num2str(Score) ])
    disp(['Final Suction Cup Grasp Score_freq = ' num2str(Score_freq) ])
    disp(['Final Suction Cup Grasp Score_amp = ' num2str(Score_amp) ])
    
                        PickMetrics(pick_point_no, 1)           = Std1;
                    PickMetrics(pick_point_no, 2:4)         = Score_metric;
                    
                    Scores_ABB      = PickMetrics(:,2);
                All_Metrics     = [ PickMetrics(:,2:4) ];
    
    
    %----------------------------------------------------------------------
    %   Overlay suction cup lip 3D point clouds in world point cloud plot
    
    figure(f1);
    tpp=text(PickPoint(1),PickPoint(2),PickPoint(3),num2str(pick_point_no));
    plot3(lipX + PickPoint(1), lipY + PickPoint(2), lipZ + PickPoint(3), 'r.');
end



%----------------------------------------------------------------------
%   Plot 3D point cloud - ALTERNATE VERSION

f3=figure;
cld_world   = plot3(worldX, worldY, worldZ, 'b');
hold on
cld_world_o = plot3(worldX', worldY', worldZ', 'r');


colormap('gray')

axis equal
xlabel('x')
ylabel('y')
title('World 3D Point Cloud')

%----------------------------------------------------------------------
%   Plot 3D point cloud in 2D Gray

f4=figure;
Zmin = min(min(world.Z));
Zmax = max(max(world.Z));

plot_Z_range = [ Zmin - 0.010 ,  Zmax + 0.010 ];

%image(rangeX, flipud(rangeY), flipud(world.Z/Zmax*64) )
image([ min(rangeX) max(rangeX) ],[ min(rangeY) max(rangeY) ], world.Z/Zmax*64 );
ax=gca;
set(ax, 'Ydir', 'normal');
colormap gray

hold on

set(f4,'DefaultLineLineWidth',3);
for pick_point_no = 1:N_pick,
    PickPoint   = all_PickPoints(pick_point_no,:);

    plot( PickPoint(1), PickPoint(2), 'y*')
end


%----------------------------------------------------------------------
%   Plot the 3D point cloud with suction cup lip point clouds REMOVED

All_Metrics_Percent     = [ (1:N_pick)'   round(All_Metrics*100) ];


S               = size(worldX);
NworldX         = S(1)*S(2);
all_Ind         = 1:NworldX;
ind_non_pick_points = setdiff( all_Ind, all_linear_Ind_Union );

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
    
    %                 t91=text(PickPoints_Bin(pick_point_no,1), ...
    %                     PickPoints_Bin(pick_point_no,2), ...
    %                     PickPoints_Bin(pick_point_no,3), ...
    %                     [ num2str(pick_point_no) ' : ' num2str(M1(pick_point_no))  ' : ' num2str(M2(pick_point_no)) ] );
    
    %             Pick_no_string  = num2str( pick_point_no );
    %             MIT_string      = num2str( round(100*All_Metrics(pick_point_no,1) ));
    %             Score_string    = num2str( round(100*All_Metrics(pick_point_no,2:5) )); % Note! ABB has std(lipZ), Score, Score_amp, Score_freq
    
    Pick_no_string  = num2str( All_Metrics_Percent( pick_point_no, 1 ) );
    Score_string    = num2str( All_Metrics_Percent( pick_point_no, 2:4 ) ); % Note! ABB has ** REMOVED std(lipZ) ***, Score, Score_amp, Score_freq
    
    text_x_offset   = Suction.RO * 0.4;
    text_y_offset   = Suction.RO * 0.3;
    text_z_offset   = 0.02;
    
    t91=text(PickPoints_Bin(pick_point_no,1) + 0, ...
        PickPoints_Bin(pick_point_no,2) + text_y_offset, ...
        PickPoints_Bin(pick_point_no,3) + text_z_offset, ...
        [ Pick_no_string ] );
    
    set(t91,'FontSize',16);
    set(t91,'FontWeight','Bold');
    set(t91,'Color','y');
    
    
    t92=text(PickPoints_Bin(pick_point_no,1) - text_x_offset, ...
        PickPoints_Bin(pick_point_no,2) - text_y_offset, ...
        PickPoints_Bin(pick_point_no,3) + text_z_offset, ...
        [  Score_string ] );
    
    set(t92,'FontSize',12);
    set(t92,'FontWeight','Bold');
    set(t92,'Color','y');
end
hold off

colormap('gray')

axis equal
xlabel('x')
ylabel('y')
title('Bin Suction Cup Lip Points : "Pick point #: ABB Score"')
xlim(plot_X_range)
ylim(plot_Y_range)
zlim(plot_Z_range)
