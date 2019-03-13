%
%   test_Objects.m
%
%   Martin Krucinski
%
%   2019-02-05 test_Wedge.m
%
%   2019-03-13 test_Objects.m
%
%   Script to test suction model with auto-generated object point clouds

in      = 0.0254;
mm      = 1e-3;
dx      = 1*mm;         % world grid step size
L       = 200*mm;       % world size
W       = 100*mm;

version = 3             % wedge version 1 - one single step
                        % 2 - two steps
                        % 3 - ramp

R       = 0.030;        % Radius
N       = 50;         % number of points along x & y directions

h       = 0.050;        % cylinder base height
A       = 0%0.0125%0.010;        % ramp height


init_Suction

slope   = Suction.max_slope_linear / (1/2);

[wedgeX, wedgeY, wedgeZ] = Wedge(R, A, N, slope, version);

%----------------------------------------------------------------------
%   Empty world

rangeX              = (-L/2):dx:(L/2);
rangeY              = (-W/2):dx:(W/2);
NX                  = length(rangeX);
NY                  = length(rangeY);

[emptyX,emptyY]     = meshgrid( rangeX , rangeY );
emptyZ             = emptyX .* 0.00;   % flat plane




%----------------------------------------------------------------------
%   Testing with wedge for now

midX        = NX/2;
midY        = NY/2;


% worldX      = wedgeX;
% worldY      = wedgeY;
% worldZ      = wedgeZ;

world.X     = worldX;
world.Y     = worldY;
world.Z     = worldZ;

wedge.X     = wedgeX;
wedge.Y     = wedgeY;
wedge.Z     = wedgeZ;



world       = addObject( empty , wedge, wedge_insert_pos );

worldY      = wedgeY;
worldZ      = wedgeZ;

%----------------------------------------------------------------------
%   Plot 3D point cloud

f1=figure;
cld_world   = plot3(worldX, worldY, worldZ, 'k.');

colormap('gray')

axis equal
xlabel('x')
ylabel('y')
title('World 3D Point Cloud')

%----------------------------------------------------------------------
%   Evaluate suction grasp score


PickPoint   = [0 0 0];  %   Wedge pick point is at origin
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

%----------------------------------------------------------------------
%   Plot 3D point cloud

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

%Score = Suction_Model_1_0( lipX, lipY, lipZ, Suction);

Score_metric = Suction_Model_1_0b( lipX, lipY, lipZ, Suction);

Score       = Score_metric(1);
Score_freq  = Score_metric(2);
Score_amp   = Score_metric(3);

if ~Lip_Complete_Flag,
    Std1 = Std_Bad;
    Score   = 0.00;
end

disp(['Final Suction Cup Grasp Score = ' num2str(Score) ])
disp(['Final Suction Cup Grasp Score_freq = ' num2str(Score_freq) ])
disp(['Final Suction Cup Grasp Score_amp = ' num2str(Score_amp) ])
