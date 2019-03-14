function wedge = Wedge2(R,h,A,slope_linear,version, dx)
%Wedge Function to generate 3D point cloud for a step wedge
%   centered around (x,y) = (0,0)
% R - cylinder radius
% h - cylinder base height
% A - step amplitude
% *** DON'T USE THIS ANY MORE, use dxN - number of point cloud points
% version = 1 - one single step
%           2 - two steps
%           3 - ramp
% dx - grid spacing

%[x,y]   = meshgrid(-R:(2*R/(N-1)):R, -R:(2*R/(N-1)):R);
[x,y]   = meshgrid(-R:dx:R, -R:dx:R);

switch version,
    case 1
        z       = h + (pi + atan2(y,x))/(2*pi)*A;
    case 2
        z       = h + (y>.0)*A;
    case 3
        ym      = A / slope_linear;
        z       = h + min(A, max(-A, (y./ym.*A)));
end

R1      = sqrt(x.^2+y.^2);
i1      = find(R1>R);
z(i1)   = zeros(size(i1));

wedge.X     = x;
wedge.Y     = y;
wedge.Z     = z;
wedge.dx    = dx;

end

