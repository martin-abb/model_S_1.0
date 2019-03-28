%
%   create_PC_from_data.m
%
%   Martin Krucinski
%
%   2019-03-28

load result_test_Objects02.mat

worldXYZ_pick_points = [ worldX_pick_points worldY_pick_points worldZ_pick_points];

pc_pick_points = pointCloud(worldXYZ_pick_points);
f1=figure;
pcshow(pc_pick_points)

worldXYZ_non_pick_points = [ worldX_non_pick_points' worldY_non_pick_points' worldZ_non_pick_points'];

pc_non_pick_points = pointCloud(worldXYZ_non_pick_points);
f2=figure;
pcshow(pc_non_pick_points)

%   add RED color to suction cup lip pick point clouds

SZ_pick_points  = size(pc_pick_points.Location);

color_pick_points   = uint8( ones(SZ_pick_points(1),1) * [255 0 0] );

pc3     = pointCloud( pc_pick_points.Location, 'color', color_pick_points );
f3=figure;
pcshow(pc3)

%   combine point clouds

%   add BLACK color to non pick point clouds

SZ_non_pick_points  = size(pc_non_pick_points.Location);

color_non_pick_points   = uint8( ones(SZ_non_pick_points(1),1) * [0 0 0] );
pc4     = pointCloud( pc_non_pick_points.Location, 'color', color_non_pick_points );

pc_all  = pcmerge( pc3 , pc4, 0.5*mm);
f4=figure;
pcshow(pc_all)

%   To be able to INDIVIDUALLY modify markers for suction cup lip point
%   clouds, plot TWO seperate point clouds in ONE plot

f5=figure;

hpc3    = pcshow(pc3);  % the pick points in red

set(hpc3.Children, 'Marker', 'o');
set(hpc3.Children, 'LineWidth', 2);

hold on
hpc4    = pcshow(pc4);  % the non pick points in black
