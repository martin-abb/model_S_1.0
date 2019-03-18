function Score_metric = Suction_Model_1_0b( lipX, lipY, lipZ, Suction)
%Suction_Model_1_0 Function implements the ABB Suction 1.0 model for
% suction cup grasping.
% Needs to be used in conjunction with function Lip_Complete_Check
%
%   Returns Score which is a normalized [0-1] metric of the probability of
%   a successful grasp at the pick point givent by the extracted suction
%   cup lip point cloud lipX, lipY, lipZ

%   Analyze lip segments in an angle increment of dtheta

debug = 0;

%*** seems I was getting NaNs in the 
% lipZ_avg_tangential
% 
% lipZ_avg_tangential =
% 
%        NaN
%    -0.0080
%        NaN
%    -0.0124
%    -0.0062
%    -0.0114
%    -0.0061
% 
% Which leads to 
% std_lipZ =
% 
%    NaN
%
%   *** SOLUTION is to either
%   1) decrease N_segments so that there are always lipZ points in that
%   angle segment
%
%   OR
%
%   try to put std = 99, avg = 0 in segments with zero points,
%   replaced the 99 later with the median of the lipZ tangential
   
N_segments      = 100;
%***N_segments      = 18;


dtheta          = 2*pi/100;

deg             = pi/180;

%*** 2019-03-18  slope_factor    = 0.2;%***1;                    % tuning factor as to how fast the score changes from 1 to 0 in the transition
slope_factor    = 0.02;%***1;                    % tuning factor as to how fast the score changes from 1 to 0 in the transition

% region from low slope (dlipZ) to
% high slope
% higher value means quicke drop-off to score 0

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
    
    if length(lipZ_segment)>0,
        lipZ_avg_tangential(i)  = mean(lipZ_segment);
        lipZ_std_tangential(i)  = std(lipZ_segment);
    else
        lipZ_avg_tangential(i)  = 0;
        lipZ_std_tangential(i)  = 99;
    end
    
    
    
end

i_99                        = find(lipZ_std_tangential == 99);
std_median                  = median(lipZ_std_tangential);
lipZ_std_tangential(i_99)   = std_median * ones(size(i_99));

%*** 2019-02-18 Added modification of lipZ_avg_tangential as well to AVOID
% SPIKES in the lipZ height and slope which can LOWER the Score_freq metric

avg_median                  = median(lipZ_avg_tangential);
lipZ_avg_tangential(i_99)   = avg_median * ones(size(i_99));

%   replicate value at index  (N_segments - 1) for the last value at
%   (N_segments)
%    **** ACTUALLY THIS IS WRONG!
%    **** REPLICATE THE 1st value instead so that we can detect a sharp
%    step right at the end!
%
%   2019-02-12  Do note that missing a step RIGHT AT THE END is UNLIKELY
%               and the above feature would make testing more difficult
%               with ramps in the middle of the profile but a step at the
%               end

%   ***** IN PROGRESS *****

lipZ_avg_tangential(N_segments)  = lipZ_avg_tangential(N_segments-1);
lipZ_std_tangential(N_segments)  = lipZ_std_tangential(N_segments-1);

%   differentiate along the tangential direction
dlipZ_avg_tangential             = [ diff(lipZ_avg_tangential) / ( dtheta * Suction.RM ) ; 0];  % LINEAR slope

%   find max slope

%***max_dlipZ           = max(abs(dlipZ_avg_tangential));   % need the max ABS slope, both negative and positive
% MK 2019-02-25 Replace with ALTERNATE metric that properly detects slope
% both BEFORE and AFTER an edge
% Per my 2019-02-25 note:
% Another option to handle this is:
% 	To handle the edge of box case better,
% 	instead of using max(dlipZ) as the max slope,
% 	use max( (dlipZ – min(dlipZ)),
% i.e. find the MAX difference in slope values…

min_dlipZ           = min(dlipZ_avg_tangential);
max_dlipZ           = max( dlipZ_avg_tangential - min_dlipZ );  % need the max slope,
                                                                % accounting for both NEG and POS slopes surrounding an edge


%***Score               = 1 - ( atan( max_dlipZ / max_slope_linear * slope_factor ) / (pi/2) );
%
%   2019-02-12  Introduce non-linear deadband, allow slopes <
%   max_slop_linear to receive a score of 100%
%
%   later, introduce the LP filtered version...
%
Score_freq              = 1 - ( atan( max( 0, max_dlipZ - Suction.max_slope_linear ) / Suction.max_slope_linear * slope_factor ) / (pi/2) );

%  2019-02-12   Include amplitude effect

std_lipZ                = std(lipZ_avg_tangential);     % overall standard deviation along tangential direction
%       Include a factor of 3 to account for better scores for larger
%       apmlitude steps with low frequency content
Score_amp               = max( 0, 1 - (std_lipZ / (3 * Suction.max_lipZ_amplitude) ) );

%   2019-02-12  Calculate combined score
Score                   = Score_amp * Score_freq;

Score_metric            = [ Score Score_freq Score_amp ];   % return these metrics

%   print for debug purposes
if debug,
    
    Score_freq
    Score_amp
    Score
    
end

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

