function [lipZSorted, lipAngles, lipZ_sort_indices ] = Sort_Lip_Z(lipX, lipY, lipZ, Suction)
%Sort_Lip_Z Sort the extracted lip Z-height values by angle for suction cup
%center


    %   Algorithm for evaluating point cloud point variations along the
    %   suction cup circle center line
    %
    %   Ideally, all points should be projected onto the circle center line,
    %   and then be evaluated based on distance along the center line,
    %   i.e.in order around the circle,distance will equal the angle
    
    %   Implementation of my suggested algorithm:
    %
    %   For all point cloud points on the suction cup lip area:
    %
    %   loop over all points, i.e. for PL_i
    %       project PL_i onto x,y plane
    %       calculate angle_i, angle along center line circle for projected
    %       point
    %   end
    %   sort all point cloud point Z-coordinates according to increasing
    %   angle_i
    %   Now, evaluate variation along the lip center line
    
    Nsize   = size(lipZ);
    
    lipAngles   = zeros(Nsize);
    angles1     = atan2(lipY, lipX);
    angles2     = unwrap(angles1);
    [Sortangles,lipZ_sort_indices]    = sort(angles2);
    lipZSorted    = lipZ(lipZ_sort_indices);        %   re-arrange height coordinates according to distance along lip centerline

end

