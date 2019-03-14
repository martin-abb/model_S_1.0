function box = Box(L, W, H, dx)
%box Function to generate 3D point cloud for a box
%   centered around (x,y) = (0,0) laying down on a flat surface, 
% L - length (x-direction)
% W - width (y-direction)
% H - box height

% dx - grid spacing

[x,y]   = meshgrid(-L/2:dx:L/2, -W/2:dx:W/2);



z       = ones(size(x))*H;

box.X     = x;
box.Y     = y;
box.Z     = z;
box.dx    = dx;

end

