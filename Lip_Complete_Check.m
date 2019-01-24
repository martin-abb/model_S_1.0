function Lip_Complete_Flag = Lip_Complete_Check( lipX, lipY, lipZ, Suction)
%Lip_Complete_Check Function to check for completeness of suction cup lip
% point cloud along the tangential direction of the lip

[CMx, CMy]  = circle(Suction.RM);

NM          = length(CMx);      % number of lip points
NSkip       = 10;               % number of points to skip during checking

Points_Found    = true;
iCM             = 1;

while Points_Found && iCM <= NM,
    Lx      = CMx(iCM);         % extract the lip point we need to check
    Ly      = CMy(iCM);
    
    DL     = sqrt((lipX - Lx).^2 + (lipY-Ly).^2);
    ind_Close_points    = find(DL <= Suction.Rcheck);
    N_Close_points      = length(ind_Close_points);
    
    if N_Close_points == 0,
        %   no close points found, hole exists in suction cup lip
        Points_Found = false;
    end
    
    iCM     = iCM + NSkip;
end

Lip_Complete_Flag = Points_Found;

end

