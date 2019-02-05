%
%   init_Suction.m
%
%   Martin Krucinski
%
%   2019-02-05
%
%   Script to initialize suction cup parameters


%----------------------------------------------------------------------
%   Define suction cup parameters
%
%in      = 0.0254;
%mm      = 1e-3;

Suction.RO      = 47/2*mm;  %20*mm %*****47/2*mm;        %   Suction cup OUTER radius
Suction.RI      = 35/2*mm;  %15*mm %***35/2*mm;        %   Suction cup INNER radius, lip width approx. 5 mm


Suction.RM      = mean([ Suction.RI Suction.RO]);    % mean suction cup lip radius
Suction.Rcheck  = (Suction.RO - Suction.RI)/2;      % radius distance to check around for lip points, half of the suction lip width
