function [worldX,worldY, worldZ] = project_D_to_XYZ(inputDepth,cameraIntrinsics,cameraPose)
%project_D_to_XYZ Project depth map into 3D point cloud

[NY,NX]=size(inputDepth)

raw_ind_x   = 1:NX;
raw_ind_y   = 1:NY;

[pixX,pixY] = meshgrid(raw_ind_x,raw_ind_y);

camX = (pixX-cameraIntrinsics(1,3)).*inputDepth/cameraIntrinsics(1,1);
camY = (pixY-cameraIntrinsics(2,3)).*inputDepth/cameraIntrinsics(2,2);
camZ = inputDepth;
worldX = cameraPose(1,1) * camX + cameraPose(1,2) * camY + cameraPose(1,3) * camZ + cameraPose(1,4);
worldY = cameraPose(2,1) * camX + cameraPose(2,2) * camY + cameraPose(2,3) * camZ + cameraPose(2,4);
worldZ = cameraPose(3,1) * camX + cameraPose(3,2) * camY + cameraPose(3,3) * camZ + cameraPose(3,4);

end

