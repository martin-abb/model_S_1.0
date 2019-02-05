%
%   plot_wedge.m
%
%   Martin Krucinski
%
%   2019-02-05
%
%   Script to generate and plot a 3D cylinder with a spiral wedge top surface

R       = 0.030;        % Radius
N       = 50;         % number of points along x & y directions

h       = 0.050;        % cylinder base height
A       = 0.010;        % ramp height

[x,y]   = meshgrid(-R:(2*R/(N-1)):R, -R:(2*R/(N-1)):R);

z       = h + (pi + atan2(y,x))/(2*pi)*A;

R1      = sqrt(x.^2+y.^2);
i1      = find(R1>R);
z(i1)   = zeros(size(i1));

f1      = figure;

surfl(x,y,z)
axis equal


