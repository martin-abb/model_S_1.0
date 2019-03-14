function cylinder = Cylinder(R,h, dx)
%Cylinder Function to generate 3D point cloud for a cylinder
%   centered around (x,y) = (0,0) laying down on a flat surface, cut in
%   half
% R - cylinder radius
% h - cylinder base height

% dx - grid spacing

[x,y]   = meshgrid(-h/2:dx:h/2, -R:dx:R);



z       = sqrt(max(0, R^2 -  y.^2));

cylinder.X     = x;
cylinder.Y     = y;
cylinder.Z     = z;
cylinder.dx    = dx;

end

