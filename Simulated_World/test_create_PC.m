%
%   test_create_PC.m
%
%   Test of creating a point cloud in ML
%
%   Martin Krucinski
%   2019-03-28
%

%   Locations
Loc2(1,1,:) = [0 0 0];
Loc2(1,2,:) = [0 0 1];
Loc2(2,1,:) = [1 0 0];
Loc2(2,2,:) = [1 1 2];
Loc2=single(Loc2);

%   Colors
Col2(1,1,:)=uint8([0 0 0]);
Col2(1,2,:)=uint8([255 0 0]);
Col2(2,1,:)=uint8([0 255 0]);
Col2(2,2,:)=uint8([0 0 255]);

pc2=pointCloud(Loc2,'Color',Col2);

f1 = figure;
hpc2=pcshow(pc2);

Flim=[-2 2];
xlim(Flim); ylim(Flim); zlim(Flim);

set(hpc2.Children, 'LineWidth', 3);
set(hpc2.Children, 'Marker', 'o');
