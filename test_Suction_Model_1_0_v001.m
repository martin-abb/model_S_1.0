%
%   test_Suction_Model_1_0_v001.m
%
%   Martin Krucinski 01/29/2019
%
%   Test script for testing the
%   Suction_Model_1_0 function

%----------------------------------------------------------------------
%   Define suction cup parameters
%
in      = 0.0254;
mm      = 1e-3;

Suction.RO      = 20*mm %*****47/2*mm;        %   Suction cup OUTER radius
Suction.RI      = 15*mm %***35/2*mm;        %   Suction cup INNER radius, lip width approx. 5 mm
Suction.RM      = mean([ Suction.RI Suction.RO]);    % mean suction cup lib radius

N_segments_test      = 1000;            % ** ?? DOES the number for testing need to be much higher than in the model evaluation ???
dtheta_test          = 2*pi/N_segments_test;

[CMx, CMy]  = circle2(Suction.RM, dtheta_test);

lipX        = CMx;
lipY        = CMy;

test_case       = 3

switch test_case,
    case 1,     % flat surface
        lipZ        = 0.011 * ones(size(lipX));
    case 2,     % step
        SL          = length(lipX);
        midSL       = round(SL/2);
        stepAmpl    = 0.005%***0.0005;
        lipZ        = [ zeros(length(1:midSL),1) ; ones(length((midSL+1):SL),1) ] * stepAmpl;
    case 3,      % smooth ramp
        
        SL          = length(lipX);
        Nmid        = 0.10 * SL;
        midSL       = SL/2;
        %   deinfe key step points
        N1          = 1;
        N2          = round(midSL - Nmid/2);
        N3          = round(midSL + Nmid/2);
        N4          = SL;
        
        stepAmpl    = 0.005%***0.0005;
        lipZ        = [ zeros((N2-N1+1),1) ; ...
                        (((N2:(N3-1))-N2)/Nmid)'* stepAmpl;
                        ones((N4-N3+2),1) * stepAmpl ];
end

Score = Suction_Model_1_0( lipX, lipY, lipZ, Suction)

2+3;