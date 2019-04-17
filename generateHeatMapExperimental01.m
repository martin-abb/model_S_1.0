%
%   generateHeatMapExperimental01.m
%
%   Martin Krucinski
%
%   2019-03-15  Inital version of separate script for heat map generation
%
%   2019-04-16  Version to run after experimental 3D PC generation from
%               RealSense camera

debug = 0;

if ~exist('tempX'),
    disp('Modifying worldX etc. from [mm] to [m].');
    tempX       = worldX;
    tempY       = worldY;
    tempZ       = worldZ;
    worldX      = tempX*mm;
    worldY      = tempY*mm;
    worldZ      = tempZ*mm;
else
    disp('Re-using worldX already converted to [m], backup in tempX etc.')
end




% Override suction cup params

% Suction.RO      = 15*mm/2
% Suction.RI      = 10*mm/2
% Suction.RM      = mean([ Suction.RO Suction.RI ])

NskipX          = 100%10;   % pixel spacing between pick points to evaluate
NskipY          = 100%10;
%   size = 3 pixels takes 334 s for the 0.6 x 0.4 m bin at dx = 1 mm,
%   27,000 points to evaluate

Score_floor     = 0.20;         % score if suction lip completely on floor, i.e z = 0.00

%----------------------------------------------------------------------
%   Generate pick points for heat map

[PickPoints,pick_set_X,pick_set_Y]    = Generate_Pick_Points(NskipX, NskipY, worldX, worldY, worldZ);

%----------------------------------------------------------------------
%
%   Process a list of pick points
%

all_linear_Ind_Union  = [];
all_linear_Ind  = {};
all_lipX        = {};
all_lipY        = {};
all_lipZ        = {};

Spick           = size(PickPoints);
N_pick          = Spick(1);
PickMetrics     = NaN(N_pick, N_metrics);   %   Matrix to store all evaluation metrics

all_Lip_Complete_Flag = zeros(N_pick,1);

disp([ 'Total number of pick points: ' num2str(N_pick) ]);


for pick_point_no = 1:N_pick,
    if mod(pick_point_no, 100) == 0,
        % home
        disp([ 'Processing pick point ' num2str(pick_point_no) '...' ]);
    end
    
    
    PickPoint   = PickPoints(pick_point_no,:);
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
    
    
    %     %----------------------------------------------------------------------
    %     %   Plot 3D point cloud for suction cup lip
    %
    %     if debug,
    %
    %         f2=figure;
    %         cld_lip   = plot3(lipX, lipY, lipZ, 'r.');
    %
    %         colormap('gray')
    %
    %         axis equal
    %         xlabel('x')
    %         ylabel('y')
    %         zlim([ -0.05 0.10 ])
    %         title('Suction Cup Lip Point Cloud')
    %
    %
    %     end
    
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
    
    
    if debug,
        disp(['Pick point number                    = ' num2str(pick_point_no) ]);
        disp(['Final Suction Cup Grasp Score        = ' num2str( round(100*Score)       ) ]);
        disp(['Final Suction Cup Grasp Score_freq   = ' num2str( round(100*Score_freq)  ) ]);
        disp(['Final Suction Cup Grasp Score_amp    = ' num2str( round(100*Score_amp)   ) ]);
        disp(' ');
    end
    

%*** this doesn't work too well with the z-height pick point specification and the lip PC being extracted
%*** from the translated PC
%***if max(max(lipZ))<0.001, % check if suction cup lip on bin floor

%**** instead, adjust with pick point height PZ
    if max(max(lipZ + PZ))<0.001, % check if suction cup lip on bin floor
        Std1            = Score_floor;
        Score_metric    = [1 1 1] * 0.20;
    end

    PickMetrics(pick_point_no, 1)           = Std1;
    PickMetrics(pick_point_no, 2:4)         = Score_metric;
    
    Scores_ABB      = PickMetrics(:,2);
    
        
    
    All_Metrics     = [ PickMetrics(:,2:4) ];
    
    
    %----------------------------------------------------------------------
    %   Overlay suction cup lip 3D point clouds in world point cloud plot
    
    if debug,
        figure(f1);
        tpp=text(PickPoint(1),PickPoint(2),PickPoint(3),num2str(pick_point_no));
        plot3(lipX + PickPoint(1), lipY + PickPoint(2), lipZ + PickPoint(3), 'r.');
    end
    
    
end


%----------------------------------------------------------------------
%   Process metrics

MX          = length(pick_set_X);
MY          = length(pick_set_Y);

M2          = PickMetrics(:,2);


Image_M2           = reshape(M2, MY, MX)*63;

%----------------------------------------------------------------------
%   Plot RAW heat maps
fh2 = figure;
image(Image_M2)
ax=gca;
set(ax, 'Ydir', 'normal');
colormap hot
title('RAW Heat Map')

