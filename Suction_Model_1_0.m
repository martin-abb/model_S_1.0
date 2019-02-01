function Score = Suction_Model_1_0( lipX, lipY, lipZ, Suction)
%Suction_Model_1_0 Function implements the ABB Suction 1.0 model for
% suction cup grasping.
% Needs to be used in conjunction with function Lip_Complete_Check
%
%   Returns Score which is a normalized [0-1] metric of the probability of
%   a successful grasp at the pick point givent by the extracted suction
%   cup lip point cloud lipX, lipY, lipZ

%   Analyze lip segments in an angle increment of dtheta

debug = 0;

N_segments      = 100;
dtheta          = 2*pi/100;

deg             = pi/180;
max_slope       = 100e-3/(45*deg); %***10e-3/(45*deg);       % max slope is 10 mm per 45 degree, *** NEEDS CALIBRATION ***
slope_factor    = 1;                    % tuning factor as to how fast the score changes from 1 to 0 in the transition
% region from low slope (dlipZ) to
% high slope

if 0,
    % **** ACTUALLY SKIP RECTANGULAR SEGMENTS IMPLEMENTATION, USE ANGLES
    % INSTEAD
    
    %   Generate points on the suction cup circles
    
    % [CMx, CMy]  = circle2(Suction.RM, dtheta);
    % [CIx, CIy]  = circle2(Suction.RI, dtheta);
    % [COx, COy]  = circle2(Suction.RO, dtheta);
    
    
    %   Approximate the suction cup lip rectangular segments
    %   by the rectangle betwen points CI and CO for that segment
    
    % lipZ_avg_tangential     = zeros(N_segments,1);      % the average
    %
    % for i=1:(N_segments-1),
    %     RIx         = CIx(i);
    %     RIx_next    = CIx(i+1);
    %     RIy         = CIy(i);
    %     RIy_next    = CIy(i+1);
    %
    %
    % [rows, cols]        = find( (R>=Suction.RI) & (R<=Suction.RO) );
    %
    
end


%   Divide suction cup lip points in segments depending on angle

angles1     = atan2(lipY, lipX);
angles2     = unwrap(angles1);
mod_angles2 = mod( angles2, 2*pi);

lipZ_avg_tangential     = zeros(N_segments,1);      % the average lipZ in the radial direction at each lip angle
lipZ_std_tangential     = zeros(N_segments,1);      % the standard deviation of lipZ in the radial direction at each lip angle

for i=1:(N_segments-1),
    theta_start = i*dtheta;
    theta_end   = (i+1)*dtheta;
    
    
    linear_indices  = find( (mod_angles2>=theta_start) & (mod_angles2<theta_end) );
    
    lipZ_segment    = lipZ( linear_indices );
    
    lipZ_avg_tangential(i)  = mean(lipZ_segment);
    lipZ_std_tangential(i)  = std(lipZ_segment);
    
    
end

%   replicate value at index  (N_segments - 1) for the last value at
%   (N_segments)
%    **** ACTUALLY THIS IS WRONG!
%    **** REPLICATE THE 1st value instead so that we can detect a sharp
%    step right at the end!


%   ***** IN PROGRESS *****

lipZ_avg_tangential(N_segments)  = lipZ_avg_tangential(N_segments-1);
lipZ_std_tangential(N_segments)  = lipZ_std_tangential(N_segments-1);

%   differentiate along the tangential direction
dlipZ_avg_tangential             = [ diff(lipZ_avg_tangential)/dtheta ; 0];

%   find max slope
max_dlipZ           = max(dlipZ_avg_tangential);

Score               = 1 - ( atan( max_dlipZ / max_slope * slope_factor ) / (pi/2) );

%--------------------------------------------------------------------------

if debug,
    figure
    set(gcf, 'DefaultLineLineWidth',3);
    plot(lipZ)
    
    figure
    set(gcf, 'DefaultLineLineWidth',3);
    plot(lipZ_avg_tangential)
    
    figure
    set(gcf, 'DefaultLineLineWidth',3);
    plot(diff(lipZ),'bo')
    
    figure
    set(gcf, 'DefaultLineLineWidth',3);
    plot(diff(lipZ_avg_tangential),'r.')
end

%--------------------------------------------------------------------------

5+3;
end

