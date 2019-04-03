%
%   Test_BlistPack_01.m
%
%   Martin Krucinski
%   martin.krucinski@us.abb.com
%
%   2019-02-21  v001    Initial version to test blister pack with
%                       ABB Suction cup model

%--------------------------------------------------------------------------
%   Read in CAD STL files

%   Blist

fname1   = 'BlistPackMarker_BicFull_v0.0.stl';
fname2   = 'LotionBottleSunScreenMR1AtScale_v0.1.stl';
fname3   = 'JuiceBox&Straw3atScale_v0.1.stl';

fname = fname3

[vBlist, fBlist, nBlist, cBlist, stltitleBlist]     = stlread(fname);



%--------------------------------------------------------------------------
%   Draw Blist

f1  = figure;

PBlist       = patch('faces', fBlist, 'vertices', vBlist);

if 0,
    set(PBlist, 'FaceColor',[0.5 0.5 0.8],'EdgeColor','k');
else
    set(PBlist,'FaceColor','b','FaceAlpha',1.0);
    set(PBlist,'EdgeColor','k','LineWidth',0.1);
end

axis equal
%zlim([-50 50])
view(3)
title(mk_str(stltitleBlist))
xlabel('x')
ylabel('y')
zlabel('z')

% cd 'Mesh_voxelisation'
% 
% %Voxelise the STL:
% [OUTPUTgrid] = VOXELISE(100,100,100,'sample.stl','xyz');
% 
% %Show the voxelised result:
% figure;
% subplot(1,3,1);
% imagesc(squeeze(sum(OUTPUTgrid,1)));
% colormap(gray(256));
% xlabel('Z-direction');
% ylabel('Y-direction');
% axis equal tight
% 
% subplot(1,3,2);
% imagesc(squeeze(sum(OUTPUTgrid,2)));
% colormap(gray(256));
% xlabel('Z-direction');
% ylabel('X-direction');
% axis equal tight
% 
% subplot(1,3,3);
% imagesc(squeeze(sum(OUTPUTgrid,3)));
% colormap(gray(256));
% xlabel('Y-direction');
% ylabel('X-direction');
% axis equal tight
% 
% cd ..

%--------------------------------------------------------------------------
%
%   Convert Rack to point cloud

addpath 'Mesh_voxelisation' -END

N       = 300;



[Blistgrid, gYBlist, gXBlist, gZBlist] = VOXELISE(N, N, N, fname);


%   Scale object to [m]
mm          = 1e-3;
if isequal(fname , 'BlistPackMarker_BicFull_v0.0.stl'),
    scale_cream = mm;
else
    scale_cream = 0.5;
end

%
% gXBlist             = gXBlist * mm;
% gYBlist             = gYBlist * mm;
% gZBlist             = gZBlist * mm;

gXBlist             = gXBlist * scale_cream;
gYBlist             = gYBlist * scale_cream;
gZBlist             = gZBlist * scale_cream;

%Blistv   = isosurface(gXBlist,gYBlist,gZBlist,Blistgrid,0.500);
f2 = figure;
isosurface(gXBlist,gYBlist,gZBlist,Blistgrid,0.500);
xlabel('x')
ylabel('y')
zlabel('z')
axis equal

%--------------------------------------------------------------------------
%
%   Plot RAW Blist point cloud

%[rows, cols]        = find( Blistgrid > 0);
%S                   = size(gXBlist);
%linearInd   = sub2ind(S, rows, cols);


%   2019-03-25 Voxelize with a mesh with the correct dimensions

dx      = 1;     % [mm]  Grid spacing
allmin  = min(vBlist);
allmax  = max(vBlist);
xmin    = allmin(1);
xmax    = allmax(1);
xR      = (xmin:dx:xmax);
NxR     = length(xR);

ymin    = allmin(2);
ymax    = allmax(2);
yR      = (ymin:dx:ymax);
NyR     = length(yR);

zmin    = allmin(3);
zmax    = allmax(3);
zR      = (zmin:dx:zmax);
NzR     = length(zR);


%[Blistgrid, gXBlist, gYBlist, gZBlist] = VOXELISE(NxR, NyR, NzR, 'BlistPackMarker_BicFull_v0.0.stl');
%[Blistgrid, gXBlist, gYBlist, gZBlist] = VOXELISE(N, N, N, 'BlistPackMarker_BicFull_v0.0.stl');

%   X & Y need to be swapped for the marker blister pack to get the correct
%   aspect ratio and dimensions

%*** ALREADY PERFORMED ABOVE
%[Blistgrid, gYBlist, gXBlist, gZBlist] = VOXELISE(N, N, N, 'BlistPackMarker_BicFull_v0.0.stl');

%--------------------------------------------------------------------------
%   Generate in XY plane



% [X, Y, Z]           = meshgrid(gXBlist, gYBlist, gZBlist);
[X, Y]              = meshgrid( gXBlist, gYBlist );
% 
% if 0,
%     IBlist               = find( Blistgrid > 0);
%     
%     BlistX                  = X(IBlist);
%     BlistY                  = Y(IBlist);
%     BlistZ                  = Z(IBlist);
%     
% else
%     
%     BlistX                  = X;
%     BlistY                  = Y;
%     BlistZ                  = Z;
%     
% end

nXB     = length(gXBlist);
nYB     = length(gYBlist);
BlistZ = zeros(nXB, nYB);
minZ    = min(min(min(Blistgrid)));

for x = 1:nXB,
    for y = 1:nYB,
        %** DOESN'T WORK BlistZ = max( Blistgrid(x, y, :) * gZBlist );
        [i,m]=find(squeeze(Blistgrid(x,y,:)) == true);
        if ~isempty(i),
            BlistZ(x,y) = gZBlist(max(i));
        else
            BlistZ(x,y) = minZ - 30*mm;
        end
    end
end

%--------------------------------------------------------------------------

f3 = figure;
%plot3(BlistX, BlistY, BlistZ, 'k.')
%***plot3(X, Y, BlistZ, 'k.')
plot3(X, Y, BlistZ, 'k.')
view(3)
xlabel('x')
ylabel('y')
zlabel('z')
axis equal

%--------------------------------------------------------------------------
%   Generate in XY plane + ROTATION along y-axis by 90 deg
%   new x = z
%   new z = -x
%   new y = y


%[Blistgrid, gYBlist, gXBlist, gZBlist] = VOXELISE(N, N, N, fname);
[Blistgrid, gZBlist, gYBlist, gXBlist] = VOXELISE(N, N, N, fname);
gZBlist = -gZBlist;


%   Scale object to [m]
mm          = 1e-3;
if isequal(fname , 'BlistPackMarker_BicFull_v0.0.stl'),
    scale_cream = mm;
else
    scale_cream = 0.5;
end

%
% gXBlist             = gXBlist * mm;
% gYBlist             = gYBlist * mm;
% gZBlist             = gZBlist * mm;

gXBlist             = gXBlist * scale_cream;
gYBlist             = gYBlist * scale_cream;
gZBlist             = gZBlist * scale_cream;


% [X, Y, Z]           = meshgrid(gXBlist, gYBlist, gZBlist);
[X, Y]              = meshgrid( gXBlist, gYBlist );
% 
% if 0,
%     IBlist               = find( Blistgrid > 0);
%     
%     BlistX                  = X(IBlist);
%     BlistY                  = Y(IBlist);
%     BlistZ                  = Z(IBlist);
%     
% else
%     
%     BlistX                  = X;
%     BlistY                  = Y;
%     BlistZ                  = Z;
%     
% end

nXB     = length(gXBlist);
nYB     = length(gYBlist);
BlistZ = zeros(nXB, nYB);
minZ    = min(min(min(Blistgrid)));

for x = 1:nXB,
    for y = 1:nYB,
        %** DOESN'T WORK BlistZ = max( Blistgrid(x, y, :) * gZBlist );
        [i,m]=find(squeeze(Blistgrid(x,y,:)) == true);
        if ~isempty(i),
            BlistZ(x,y) = gZBlist(max(i));
        else
            BlistZ(x,y) = minZ - 30*mm;
        end
    end
end

%--------------------------------------------------------------------------

f3 = figure;
%plot3(BlistX, BlistY, BlistZ, 'k.')
%***plot3(X, Y, BlistZ, 'k.')
plot3(X, Y, BlistZ, 'k.')
view(3)
xlabel('x')
ylabel('y')
zlabel('z')
axis equal


%--------------------------------------------------------------------------

f4 = figure;
%plot3(BlistX, BlistY, BlistZ, 'k.')
%***plot3(X, Y, BlistZ, 'k.')
plot3(Y, Z, BlistX, 'k.')
view(3)
xlabel('y')
ylabel('z')
zlabel('x')
axis equal


%--------------------------------------------------------------------------
%   Generate in XZ plane

% [X, Y, Z]           = meshgrid(gXBlist, gYBlist, gZBlist);
[X, Z]              = meshgrid( gXBlist, gZBlist );
% 
% if 0,
%     IBlist               = find( Blistgrid > 0);
%     
%     BlistX                  = X(IBlist);
%     BlistY                  = Y(IBlist);
%     BlistZ                  = Z(IBlist);
%     
% else
%     
%     BlistX                  = X;
%     BlistY                  = Y;
%     BlistZ                  = Z;
%     
% end


nXB     = length(gXBlist);
nZB     = length(gZBlist);
BlistY = zeros(nXB, nZB);
minY    = min(min(min(Blistgrid)));

for x = 1:nXB,
    for z = 1:nZB,
        %** DOESN'T WORK BlistZ = max( Blistgrid(x, y, :) * gZBlist );
        [i,m]=find(squeeze(Blistgrid(x,:,z)) == true);
        if ~isempty(i),
            % normall the top surface 
            % BlistY(x,z) = gYBlist(max(i));
            % for juice box, grab bottom surface
            BlistY(x,z) = gYBlist(min(i));
        else
            BlistY(x,z) = minY - 30*mm;
        end
    end
end

%--------------------------------------------------------------------------

f5 = figure;
%plot3(BlistX, BlistY, BlistZ, 'k.')
%***plot3(X, Y, BlistZ, 'k.')
plot3(X, Z, BlistY, 'k.')
view(3)
xlabel('x')
ylabel('z')
zlabel('y')
axis equal


