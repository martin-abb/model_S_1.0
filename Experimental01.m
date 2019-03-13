%
%   Experimental01.m
%
%   Using Suction 1.0 model to compare MIT pick points with ABB Suction 1.0 model
%   Martin Krucinski
%
%
%   2019-02-28      Comparison05.m Including information of WHICH MIT pick
%                   point was used AND if it was successful

%   2019-03-13      Experimental01.m
%                   Perform analysis on images from an experimental
%                   RealSense camera

%   Define running mode

script_flag = 1;
plot_PC_overlay = 0;


%   Set up looping over directories


base_dir_local  = 'C:\Users\USMAKRU\OneDrive - ABB\2018\Logistics\Suction Cup Modeling\Suction_1.0_ML';
data_dir_local  = 'C:\Users\USMAKRU\Documents\ABB Local\Intel RealSense\Tests_2019_03_13';
base_dir        = base_dir_local;
data_dir        = data_dir_local;

cd(data_dir)

%----------------------------------------------------------------------
%   Define suction cup parameters

in      = 0.0254;
mm      = 1e-3;

init_Suction

%----------------------------------------------------------------------
%   Define bin parameters

init_Bin

%----------------------------------------------------------------------
%   Define grasp quality evaluation parameters
%
Std_Bad      = 0.150; %***0.030;           % standard deviation for pick points to be flagged as BAD

%----------------------------------------------------------------------
%   Load sensor data images for testing
%

cameraIntrinsics    = load('test-camera-intrinsics.txt','ascii');
% cameraPose          = load('test-camera-pose.txt','ascii');
% 0.9999    0.0115    0.0050    0.1189
% 0.0062   -0.8006    0.5991   -0.1604
% 0.0109   -0.5990   -0.8006    0.7012
% 0         0         0    1.0000

cameraPos           = [ 0.30; 0.03 ; 1.02];     % test camera position
%cameraR             = [ 1  0  0 ; 0  1  0 ; 0  0  -1 ]; % test camera rotation, pointing downwards
cameraR             = [ 1  0  0 ; 0  -1  0 ; 0  0  -1 ]; % test camera rotation, pointing downwards
%cameraR             = [ 1  0  0 ; 0  -1  0 ; 0  0  -11 ]; % test camera rotation, pointing downwards
cameraPose          = [ cameraR   cameraPos; 0 0 0 1];


depth_scale     = 1/10^4;

N_metrics           = 4;    % number of metrics to record for each pick point

img_name        = '0004_';

img_raw_color   = imread([ data_dir '\' img_name 'Color.png' ]);
img_raw_depth   = imread([ data_dir '\' img_name 'Depth.png' ]);

%----------------------------------------------------------------------
%   Convert depth image


%inputDepth      = double(img_raw_depth)*depth_scale;
x_min           = 0.8;
x_max           = 1.2;

inputDepth      = double(img_raw_depth(:,:,1))/255 * (x_max - x_min) + x_min;

%----------------------------------------------------------------------
%  Plot the loaded image data




%----------------------------------------------------------------------
f1=figure;
image(img_raw_color);
title('Raw Color Image')
axis equal

%----------------------------------------------------------------------
f2=figure;
image(img_raw_depth);
colormap('gray')
title('Raw Depth Image')
axis equal


%----------------------------------------------------------------------
%
%   Filter depth image to move glitch points at d = 0.0 further back, to
%   e.g. 2.0 m
%
%   Filter the inputDepth image, not the img_raw_depth



d_move_back     = 1.4;%8.0;   %   Move back distance
inputDepth_filtered     = filter_depth_exp( inputDepth, d_move_back);

% img_filtered_depth = inputDepth_filtered / depth_scale;      % create filtered depth image map
% 
% if 0,
%     %----------------------------------------------------------------------
%     f3=figure;
%     image(img_filtered_depth/256);
%     colormap('gray')
%     title('Filtered Depth Image (glitch filtered)')
% end


%----------------------------------------------------------------------

[worldX,worldY, worldZ] = project_D_to_XYZ(inputDepth_filtered,cameraIntrinsics,cameraPose);


%----------------------------------------------------------------------

if 1,
    
    f5=figure;
    cld_world   = plot3(worldX, worldY, worldZ, 'k.');
    
    colormap('gray')
    
    axis equal
    xlabel('x')
    ylabel('y')
    title('World 3D Point Cloud')
    
end

cd(base_dir)

PickPoints_Bin_RAW = tote_file_contents(:,1:3);

PickPoints   = PickPoints_Bin_RAW(:, 1:3)*mm;


PickPoints_Bin          = PickPoints;

plot_X_range = [ Bin.X - 0.1 Bin.X + Bin.BW + 0.1 ];
plot_Y_range = [ Bin.Y - 0.1 Bin.Y + Bin.BL + 0.1 ];
plot_Z_range = [ Bin.Z - 0.1 Bin.Z + Bin.BH + 0.1 ];

%----------------------------------------------------------------------

if 0,
    
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

%***    Lip_Complete_Flag = true;
Lip_Complete_Flag = Lip_Complete_Check( lipX, lipY, lipZ, Suction);

all_Lip_Complete_Flag(pick_point_no) = Lip_Complete_Flag;

% Model_1_0b now returns    Score_metric            = [ Score Score_freq Score_amp ];
Score_metric = Suction_Model_1_0b( lipX, lipY, lipZ, Suction);
Score        = Score_metric(1);

if ~Lip_Complete_Flag,
    Std1 = Std_Bad;
    Score   = 0.00;
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
        Std1    = Std_Bad;
        Score   = 0.00;
    end
else
    
    %   PickPoint outside bin
    
    if ~( (PX >= Bin.X_min) & (PX <= Bin.X_max) & ...
            (PY >= Bin.Y_min) & (PY <= Bin.Y_max) & ...
            (PZ >= Bin.Z_min) & (PZ <= Bin.Z_max) ),
        Std1    = Std_Bad;
        Score   = 0.00;
    end
    
end



%----------------------------------------------------------------------
%   Sort lipZ values along lip angle

[lipZSorted, lipAngles, lipZ_sort_indices ] = Sort_Lip_Z(lipX, lipY, lipZ, Suction);


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
if ~script_flag || plot_select_flag,
    f8      = figure;
    
    Score_metric_Percent     = [ round(Score_metric*100) ];
    subplot(221)
    plot3(lipX + PickPoint(1), lipY + PickPoint(2), lipZ + PickPoint(3), 'r.');
    t81=text(PickPoint(1),PickPoint(2),PickPoint(3),num2str(pick_point_no));
    xlabel('x')
    ylabel('y')
    view(3)
    title([ 'Pick Point # ' num2str(pick_point_no) ' ABB [Score, Freq, Amp] = ' num2str(Score_metric_Percent) ])
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

if 0,
    if plot_PC_overlay,
        figure(f5);
        hold on
        plot3(lipX + PickPoint(1), lipY + PickPoint(2), lipZ + PickPoint(3), 'r.');
        t51=text(PickPoint(1),PickPoint(2),PickPoint(3),num2str(pick_point_no));
    end
end


%----------------------------------------------------------------------
%   Save metrics for this pick point

PickMetrics(pick_point_no, 1)           = Std1;
PickMetrics(pick_point_no, 2:4)         = Score_metric;


% for DEBUG purposes
if pick_point_no==8,
    2+3;
end



%----------------------------------------------------------------------
%   Process metrics
%
%   Combine MIT & ABB metrics

Scores_ABB      = PickMetrics(:,2);
All_Metrics     = [ Scores_MIT  PickMetrics(:,2:4) ];

Success_column      = (NA_number) * ones(N_pick,1);    % Here NA_number means NOT APPLICABLE / NOT SELECTED
Error_column        = (NA_number) * ones(N_pick,1);    % Here NA_number means NOT APPLICABLE / NOT SELECTED
if Pick_Success,        % fill in a value of TRUE in the SELECTED Pick point ROW if pick was succesful
    Success_column(PickPoint_index)    = Pick_Success;
end
if Error_code_index>0,
    Error_column(PickPoint_index) = Error_code;
    %   Remove the NA_number formatting number for a Failed
    %   pick, replace with 0
    if Error_code ~= 0,
        Success_column(PickPoint_index)    = 0;
    end
    
end


All_Metrics_Percent     = [ (1:N_pick)'   round(All_Metrics*100) Success_column Error_column ];



disp('Metrics: MIT and ABB Suction 1.0')
%   Layout on screen below, match indentations for column
%   headers
%          #   MIT Score  Freq   Amp   Success
%                3          84          64          66          97       -9999
disp('           #         MIT         ABB        Freq         Amp     Success       Error')
disp(All_Metrics_Percent)

%   Re-sort metrics based on ABB score

[Sorted_ABB_Score, i_sort_ABB]  = sort(Scores_ABB,'descend');

All_Metrics_Percent_ABB_sorted = All_Metrics_Percent(i_sort_ABB, :);

disp('Metrics: MIT and ABB Suction 1.0 sorted by ABB score')
disp('           #         MIT         ABB        Freq         Amp     Success       Error')
disp(All_Metrics_Percent_ABB_sorted)



%***    if plot_PC_overlay && ~script_flag,
if plot_PC_overlay,
    
    
    %----------------------------------------------------------------------
    %   Prepare the 3D point cloud with suction cup lip point clouds REMOVED
    
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
        MIT_string      = num2str( All_Metrics_Percent( pick_point_no, 2 ) );
        Score_string    = num2str( All_Metrics_Percent( pick_point_no, 3:5 ) ); % Note! ABB has ** REMOVED std(lipZ) ***, Score, Score_amp, Score_freq
        
        text_x_offset   = Suction.RO * 0.4;
        text_y_offset   = Suction.RO * 0.3;
        text_z_offset   = 0.02;
        
        t91=text(PickPoints_Bin(pick_point_no,1) + 0, ...
            PickPoints_Bin(pick_point_no,2) + text_y_offset, ...
            PickPoints_Bin(pick_point_no,3) + text_z_offset, ...
            [ Pick_no_string ] );
        
        set(t91,'FontSize',16);
        set(t91,'FontWeight','Bold');
        set(t91,'Color','b');
        
        
        t92=text(PickPoints_Bin(pick_point_no,1) - text_x_offset, ...
            PickPoints_Bin(pick_point_no,2) - text_y_offset, ...
            PickPoints_Bin(pick_point_no,3) + text_z_offset, ...
            [  MIT_string ' : ' Score_string ] );
        
        set(t92,'FontSize',12);
        set(t92,'FontWeight','Bold');
        set(t92,'Color','b');
    end
    hold off
    
    colormap('gray')
    
    axis equal
    xlabel('x')
    ylabel('y')
    title('Bin Suction Cup Lip Points : "Pick point #: MIT Score : ABB Score"')
    xlim(plot_X_range)
    ylim(plot_Y_range)
    zlim(plot_Z_range)
    
end

%----------------------------------------------------------------------
%   Add MIT & ABB Scores to RGB image with pick points
%----------------------------------------------------------------------


figure(f4);
axis equal

%***    title('MIT Princeton Pick Points with "MIT Score : ABB Score" (RED - disagreement)')
title(mk_str([ totename '    MIT : ABB [Score, Freq, Amp] (RED - delta > 30%)' ]))
for pick_point_no = 1:N_pick,
    
    MIT_Score       = round(100*All_Metrics(pick_point_no,1) );
    ABB_Score       = round(100*All_Metrics(pick_point_no,3) );
    
    Pick_no_string  = num2str( All_Metrics_Percent( pick_point_no, 1 ) );
    MIT_string      = num2str( All_Metrics_Percent( pick_point_no, 2 ) );
    ABB_string    = num2str( All_Metrics_Percent( pick_point_no, 3:5 ) );
    % Note! ABB has Score_metric            = [ Score Score_freq Score_amp ];
    
    
    Score_threshold     = 30;
    if abs(MIT_Score - ABB_Score)>= Score_threshold,
        score_color = 'r';
    else
        score_color = 'b';
    end
    
    text_x_offset   = 20;
    text_y_offset   = 0;
    
    t92=text(pixel_file_contents(pick_point_no,1) + text_x_offset, ...
        pixel_file_contents(pick_point_no,2) - text_y_offset, ...
        [  MIT_string ' : ' ABB_string ] );
    
    set(t92,'FontSize',10);
    set(t92,'FontWeight','Bold');
    set(t92,'Color',score_color);
end


%***    saveas(f4, [ logdir '\' img_name 'compare.png']);
img_file_name   = [ img_name 'compare.png'];
if exist(img_file_name)==2,
    delete(img_file_name);
end
saveas(f4,img_file_name );


%   prepare for processing of next folder
cd(base_dir)

