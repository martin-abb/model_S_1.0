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

Suction.RO      = 10*mm%47/2*mm;  %20*mm %*****47/2*mm;        %   Suction cup OUTER radius
Suction.RI      = 7*mm%35/2*mm;  %15*mm %***35/2*mm;        %   Suction cup INNER radius, lip width approx. 5 mm


Suction.RM      = mean([ Suction.RI Suction.RO]);    % mean suction cup lip radius
Suction.Rcheck  = (Suction.RO - Suction.RI)/2;      % radius distance to check around for lip points, half of the suction lip width

%-------------------------------------------------------------------------
%
%   Max slope parameter
%
%
%max_slope       = 100e-3/(45*deg); %***10e-3/(45*deg);       % max slope is 10 mm per 45 degree, *** NEEDS CALIBRATION ***
Suction.max_slope_linear = 0.55;        % from calibration experiments on 02/01/2019, max linear slope

Suction.max_lipZ_amplitude  = Suction.RO * 2;       % max allowable Z-axis amplitude for which Score decreases to 0
