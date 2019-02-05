function [x,y,z] = Wedge(R,A,N,version)
%Wedge Function to generate 3D point cloud for a step wedge
%   centered around (x,y) = (0,0)
% R - cylinder radius
% A - step amplitude
% N - number of point cloude points
% version = 1 - one single step
%           2 - two steps
%   cylinder nominal height is 50 mm
h       = 0.050;        % cylinder base height
[x,y]   = meshgrid(-R:(2*R/(N-1)):R, -R:(2*R/(N-1)):R);

switch version,
    case 1
        z       = h + (pi + atan2(y,x))/(2*pi)*A;
    case 2
        z       = h + (y>.0)*A;
end

R1      = sqrt(x.^2+y.^2);
i1      = find(R1>R);
z(i1)   = zeros(size(i1));


end

