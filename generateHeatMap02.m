%
%   generateHeatMap02.m
%
%   Martin Krucinski
%
%   2019-03-15  Inital version of separate script for heat map generation
%   2019-03-19  v02 - utilizes Suction_Model_1_0c to include mass & torque
%   effects

load result_test_Objects03

% Override suction cup params

% Suction.RO      = 15*mm/2
% Suction.RI      = 10*mm/2
% Suction.RM      = mean([ Suction.RO Suction.RI ])

NskipX          = 3%10;   % pixel spacing between pick points to evaluate
NskipY          = 3%10;
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

%   set up empty object to use everywhere where there are NO objects in
%   order for suction cup model to evaluate

empty_object.L      = 0.150;
empty_object.W      = 0.050;
empty_object.H      = 0.040;
empty_object        = Box( empty_object.L , empty_object.W , empty_object.H , dx);
% re-insert over-written fields, re-factor this later
empty_object.L      = 0.150;
empty_object.W      = 0.050;
empty_object.H      = 0.040;

empty_object.m        = 0.000;                  % [kg]  object mass
empty_object.w        = empty_object.m * g;     % [N]   object weight
empty_object.CM_x     = 0.000;                  % [m]   relative CM position
empty_object.CM_y     = 0.000;
empty_object.CM_z     = empty_object.H/2;
empty_object.CM       = [ empty_object.CM_x  empty_object.CM_y  empty_object.CM_z  ]';

empty_object.insert_pos     = [0.000 ; 0.000 ; 0.000];



for pick_point_no = 1:N_pick,
    if mod(pick_point_no, 100) == 0,
        % home
        disp([ 'Processing pick point ' num2str(pick_point_no) '...' ]);
    end
    
    
    PickPoint   = PickPoints(pick_point_no,:)';
    PickPoint_Bin   = PickPoint;
    
    %----------------------------------------------------------------------
    %   Evaluate suction grasp score
    
    
    
    
    PX = PickPoint(1);
    PY = PickPoint(2);
    PZ = PickPoint(3);
    
    %----------------------------------------------------------------------
    %   Determine WHICH object is being picked up
    
    if ( PX>=0.100 && PX<=0.250 && PY>=0.050 && PY<=0.100 )
        Object      = box1;%all_PickObjects{1};
    elseif ( PX>=0.100 && PX<=0.250 && PY>=0.150 && PY<=0.200 )
        Object      = box2;%all_PickObjects{2};
    elseif ( PX>=0.100 && PX<=0.250 && PY>=0.250 && PY<=0.300 )
        Object      = box3;%all_PickObjects{3};
    elseif ( PX>=0.350 && PX<=0.500 && PY>=0.050 && PY<=0.100 )
        Object      = box4;%all_PickObjects{4};
    elseif ( PX>=0.350 && PX<=0.500 && PY>=0.150 && PY<=0.200 )
        Object      = box5;%all_PickObjects{5};
    elseif ( PX>=0.350 && PX<=0.500 && PY>=0.250 && PY<=0.300 )
        Object      = box6;%all_PickObjects{6};
    else
        Object = empty_object;
    end
    
        
    %             %       Box with CM in the middle
    %
    %         box1.L      = 0.150;
    %         box1.W      = 0.050;
    %         box1.H      = 0.040;
    %         box1        = Box( box1.L , box1.W , box1.H , dx);
    %         % re-insert over-written fields, re-factor this later
    %         box1.L      = 0.150;
    %         box1.W      = 0.050;
    %         box1.H      = 0.040;
    %         box1_insert_pos     = [0.100 ; 0.050 ; 0];
    %
    %         box1.m        = 0.300;                  % [kg]  object mass
    %         box1.w        = box1.m * g;             % [N]   object weight
    %         box1.CM_x     = 0.000;                  % [m]   relative CM position
    %         box1.CM_y     = 0.000;
    %         box1.CM_z     = box1.H/2;
    %         box1.CM       = [ box1.CM_x  box1.CM_y  box1.CM_z  ];
    %
    %         %       Box with CM one quarter of the length to the left
    %
    %         box2          = box1;
    %         box2.CM_x     = -box2.L/4;                  % [m]   relative CM position
    %         box2.CM       = [ box2.CM_x  box2.CM_y  box2.CM_z  ];
    %         box2_insert_pos     = [0.100 ; 0.150 ; 0];
    %
    %         %       Box with CM at the left edge
    %
    %         box3          = box1;
    %         box3.CM_x     = -box3.L/2;                  % [m]   relative CM position
    %         box3.CM       = [ box3.CM_x  box3.CM_y  box3.CM_z  ];
    %         box3_insert_pos     = [0.100 ; 0.250 ; 0];
    %
    %
    
    %----------------------------------------------------------------------
    % Set the pick pose
    
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
    
    
    % From test_Objects03
%         Object      = all_PickObjects{pick_point_no};
%     
%     ObjectMidPoint  = Object.insert_pos + [ Object.L/2 ; Object.W/2 ; Object.H/2 ];
%     
%     
%     Score_metric = Suction_Model_1_0c( lipX, lipY, lipZ, Suction, Object, PickPoint - ObjectMidPoint);
%     
    ObjectMidPoint  = Object.insert_pos + [ Object.L/2 ; Object.W/2 ; Object.H/2 ];
    
    Score_metric = Suction_Model_1_0c( lipX, lipY, lipZ, Suction, Object, PickPoint - ObjectMidPoint );
    
    Score       = Score_metric(1);
    Score_freq  = Score_metric(2);
    Score_amp   = Score_metric(3);
    
    Std1=std(lipZ);        % Metric 1 - total Z-height standard deviation
    
    if ~Lip_Complete_Flag,
        Std1 = Std_Bad;
        Score   = 0.00;
    end
    
%     
%     if debug,
%         disp(['Pick point number                    = ' num2str(pick_point_no) ]);
%         disp(['Final Suction Cup Grasp Score        = ' num2str( round(100*Score)       ) ]);
%         disp(['Final Suction Cup Grasp Score_freq   = ' num2str( round(100*Score_freq)  ) ]);
%         disp(['Final Suction Cup Grasp Score_amp    = ' num2str( round(100*Score_amp)   ) ]);
%         disp(' ');
%     end
    

%*** this doesn't work too well with the z-height pick point specification and the lip PC being extracted
%*** from the translated PC
%***if max(max(lipZ))<0.001, % check if suction cup lip on bin floor

%**** instead, adjust with pick point height PZ
    if max(max(lipZ + PZ))<0.001, % check if suction cup lip on bin floor
        Std1            = 0.20;
        Score_metric    = [1 1 1 1 1] * 0.20;
    end

    PickMetrics(pick_point_no, 1)           = Std1;
    PickMetrics(pick_point_no, 2:6)         = Score_metric;
    
    Scores_ABB      = PickMetrics(:,2);
    
        
    
    All_Metrics     = [ PickMetrics(:,2:6) ];
    
    
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
Image_M2_Large      = Image_Block_Expand(Image_M2, NskipX, NskipY);

%----------------------------------------------------------------------
%   Plot RAW heat maps
fh2 = figure;
image( rangeX, rangeY, Image_M2_Large)
axis equal
ax=gca;
set(ax, 'Ydir', 'normal');
colormap hot
title('RAW Heat Map')

%---------------------------
%   Plot heat map with CENTER OF MASS symbols

fh3 = figure;
Image_M2_Large_CM   = Image_M2_Large;

for o = 1:length(all_PickObjects),
    Object     = all_PickObjects{o};
    PCM_rel    = Object.CM;       % relative CM position  in [m]
    ObjectMidPoint  = Object.insert_pos + [ Object.L/2 ; Object.W/2 ; Object.H/2 ];
    PCM     = PCM_rel + ObjectMidPoint;
    
    PX      = PCM(1);
    PY      = PCM(2);
    
    % convert position to pixels
    %*** this code for SMALL not LARGE heat map    ICMX    = floor( (PX - min(rangeX)) / dx / NskipX ) + 1;
    %***    ICMY_temp    = floor( (PY - min(rangeY)) / dx / NskipY ) + 1;
    
    %   Code for LARGE heat map
    ICMX    = floor( (PX - min(rangeX)) / dx ) + 1;
%     ICMY_temp    = floor( (PY - min(rangeY)) / dx ) + 1;
%     ICMY    = NY - ICMY_temp + 1;   % flip Y-axis
    %   no Y-axis flipping required!
    ICMY    = floor( (PY - min(rangeY)) / dx ) + 1;

    %**    ICM     = [ ICMX ICMY ];
    ICM     = [ ICMX ICMY ];  
    Image_M2_Large_CM = draw_CM(Image_M2_Large_CM, ICM,10);
    
end

image( rangeX, rangeY, Image_M2_Large_CM)
axis equal
ax=gca;
set(ax, 'Ydir', 'normal');
colormap hot
title('RAW Heat Map with CM locations')
