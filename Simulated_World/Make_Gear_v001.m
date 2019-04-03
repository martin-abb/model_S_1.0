%
%   Make_Gear_v001.m
%
%   Martin Krucinski
%   martin.krucinski@us.abb.com
%
%   2019-02-21  v001    Initial version to test feasiblity of Matlab gear
%                       generation

%--------------------------------------------------------------------------
%   Read in CAD STL files

%   Rack
[vRack, fRack, nRack, cRack, stltitleRack]     = stlread('Rack.STL');

%   Cylinder
[vCylinder, fCylinder, nCylinder, cCylinder, stltitleCylinder]     = stlread('Cylinder.STL');

%--------------------------------------------------------------------------
%   Draw Rack

f1  = figure;

PRack       = patch('faces', fRack, 'vertices', vRack);

if 0,
    set(PRack, 'FaceColor',[0.5 0.5 0.8],'EdgeColor','k');
else
    set(PRack,'FaceColor','b','FaceAlpha',1.0);
    set(PRack,'EdgeColor','k','LineWidth',0.1);
end

axis equal
%zlim([-50 50])
view(3)
title(mk_str(stltitleRack))


%--------------------------------------------------------------------------
%   Draw Cylinder

f2  = figure;

PCylinder       = patch('faces', fCylinder, 'vertices', vCylinder);
set(PCylinder, 'FaceColor',[0.5 0.5 0.8],'EdgeColor','k');
axis equal
%zlim([-100 100])
view(3)
title(mk_str(stltitleCylinder))

%--------------------------------------------------------------------------

cd 'Mesh_voxelisation'

%Voxelise the STL:
[OUTPUTgrid] = VOXELISE(100,100,100,'sample.stl','xyz');

%Show the voxelised result:
figure;
subplot(1,3,1);
imagesc(squeeze(sum(OUTPUTgrid,1)));
colormap(gray(256));
xlabel('Z-direction');
ylabel('Y-direction');
axis equal tight

subplot(1,3,2);
imagesc(squeeze(sum(OUTPUTgrid,2)));
colormap(gray(256));
xlabel('Z-direction');
ylabel('X-direction');
axis equal tight

subplot(1,3,3);
imagesc(squeeze(sum(OUTPUTgrid,3)));
colormap(gray(256));
xlabel('Y-direction');
ylabel('X-direction');
axis equal tight

cd ..

%--------------------------------------------------------------------------
%
%   Convert Rack to point cloud

addpath 'Mesh_voxelisation' -END

N       = 100;
[Rackgrid, gXRack, gYRack, gZRack] = VOXELISE(N, N, N, 'Rack.STL');


%Rackv   = isosurface(gXRack,gYRack,gZRack,Rackgrid,0.500);
isosurface(gXRack,gYRack,gZRack,Rackgrid,0.500);

%--------------------------------------------------------------------------
%
%   Plot RAW Rack point cloud

%[rows, cols]        = find( Rackgrid > 0);
%S                   = size(gXRack);
%linearInd   = sub2ind(S, rows, cols);

IRack               = find( Rackgrid > 0);
[X, Y, Z]           = meshgrid(gXRack, gYRack, gZRack);

RackX                  = X(IRack);
RackY                  = Y(IRack);
RackZ                  = Z(IRack);

figure
plot3(RackX, RackY, RackZ, 'k.')
view(90,0)
xlabel('x')
ylabel('y')
zlabel('z')


%--------------------------------------------------------------------------
%
%   Convert cylinder

[Cylindergrid, gXCylinder, gYCylinder, gZCylinder] = VOXELISE(N, N, N, 'Cylinder.STL');

ICylinder               = find( Cylindergrid > 0);
[X, Y, Z]               = meshgrid(gXCylinder, gYCylinder, gZCylinder);

CylinderX                  = X(ICylinder);
CylinderY                  = Y(ICylinder);
CylinderZ                  = Z(ICylinder);

figure
plot3(CylinderX, CylinderY, CylinderZ, 'b.')
view(90,0)
xlabel('x')
ylabel('y')
zlabel('z')



%--------------------------------------------------------------------------
%   Plot both Rack and Cylinder point clouds in the same plot
figure
plot3(RackX, RackY, RackZ, 'k.')
hold on
plot3(CylinderX, CylinderY, CylinderZ, 'b.')

view(3)

xlabel('x')
ylabel('y')
zlabel('z')


