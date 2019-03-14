function sphere = Sphere(R, dx)
% Sphere Function to generate 3D point cloud for a half-sphere
%   centered around (x,y) = (0,0) laying down on a flat surface, cut in
%   half
% R - sphere radius

% dx - grid spacing

[x,y]   = meshgrid(-R:dx:R, -R:dx:R);



z       = sqrt(max(0, R^2 - x.^2 - y.^2));

sphere.X     = x;
sphere.Y     = y;
sphere.Z     = z;
sphere.dx    = dx;

end

