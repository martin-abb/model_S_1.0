%
%   init_Bin.m
%
%   Martin Krucinski
%
%   2019-02-15  Initialize bin parameters, needs to run after
%   init_Suction.m

Bin.BW          = 375*mm;
Bin.BL          = 550*mm;
Bin.BH          = 220*mm;       % Bin dimensions Width x Length, Height, Origin
Bin.X           = 0.000;        % Bin origin
Bin.Y           = 0.000;
Bin.Z           = 0.000;

% Bin.dx          = 0.1;          % tolerance around bin limits to use for pick point filtering
% Bin.dy          = 0.1;
% Bin.dz          = 0.1;

Bin.dx          = Suction.RO;          % tolerance around bin limits to use for pick point filtering
Bin.dy          = Suction.RO;
%***Bin.dz          = 0.020;%****0.01;
% 2019-02-15 Increase Z-height to be able to handle very full bins that
% have pick points with Z-coordinates of 205, 228, 230 mm!!
Bin.dz          = -0.100;%****0.01;

%   The ranges below were used for point cloud axis limits, extend BEYOND
%   bin
% Bin.X_min                 = Bin.X - Bin.dx;
% Bin.X_max                 = Bin.X + Bin.BW + Bin.dx;
%
% Bin.Y_min                 = Bin.Y - Bin.dy;
% Bin.Y_max                 = Bin.Y + Bin.BL + Bin.dy;
%
% Bin.Z_min                 = Bin.Z - Bin.dz;
% Bin.Z_max                 = Bin.Z + Bin.BH + Bin.dz;

%   Use these ranges for checking if pick point inside bin, REDUCE the bin
%   size instead of extending

Bin.X_min                 = Bin.X + Bin.dx;
Bin.X_max                 = Bin.X + Bin.BW - Bin.dx;

Bin.Y_min                 = Bin.Y + Bin.dy;
Bin.Y_max                 = Bin.Y + Bin.BL - Bin.dy;

Bin.Z_min                 = Bin.Z + Bin.dz;
Bin.Z_max                 = Bin.Z + Bin.BH - Bin.dz;