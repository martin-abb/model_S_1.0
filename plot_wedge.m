%
%   plot_wedge.m
%
%   Martin Krucinski
%
%   2019-02-05
%
%   Script to generate and plot a 3D cylinder with a spiral wedge top surface

version = 1;            % 1 - one single step
                        % 2 - two steps

R       = 0.030;        % Radius
N       = 50;         % number of points along x & y directions

h       = 0.050;        % cylinder base height
A       = 0.010;        % ramp height

[x,y,z]     = Wedge(R, A, N, version);

f1      = figure;

surfl(x,y,z)
axis equal


