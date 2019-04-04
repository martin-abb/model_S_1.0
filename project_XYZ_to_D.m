function [pixelX, pixelY, Depth] = project_XYZ_to_D(worldX,worldY, worldZ,cameraIntrinsics,cameraPose, imageSize, flip)
%project_D_to_XYZ Project depth map into 3D point cloud
%
%   flip = true     : image y-coordinate with a upper left corner origin and y-axis pointing down
%   filp = false    : image y-coordinate with a lower left corner origin and y-axis pointing up

[NY,NX]     = size(worldX);
Iheight     = imageSize(2);     % image size is [ width , height ]

CIM         = cameraIntrinsics;
CEM         = cameraPose;
C           = CIM * CEM ;

raw_ind_x   = 1:NX;
raw_ind_y   = 1:NY;

for x = raw_ind_x,
    for y = raw_ind_y,
        W       = [ worldX(y,x) ; worldY(y,x) ; worldZ(y,x) ];
        Wext    = [ W ; 1 ];
        
        ptilde  = C*Wext;
        
        Depth   = ptilde(3);
        %        p       = round( ptilde / Depth );
        p       = ( ptilde / Depth );
        
        pixelX          = p(1);
        pixelY_temp     = p(2);    % y-coordinate with a lower left corner origin and y-axis pointing up
        
        
        %disp('Pixel coordinates')
        %disp([ 'pixx  = ' num2str(pixx) ', pixy = ' num2str(pixy) ', Z = ' num2str(Z) ])
        
        if flip,
            pixelY   = Iheight - pixelY_temp;   % y-coordinate with a upper left corner origin and y-axis pointing down
        else
            pixelY   = pixelY_temp;   % y-coordinate with a upper left corner origin and y-axis pointing down
        end
    end
end

