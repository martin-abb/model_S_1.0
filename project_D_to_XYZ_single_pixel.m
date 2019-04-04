function [worldX,worldY, worldZ] = project_D_to_XYZ_single_pixel(px, py, inputDepth,cameraIntrinsics,cameraPose)
%project_D_to_XYZ Project single depth point into 3D point

pixX    = px;
pixY    = py;

camX = (pixX-cameraIntrinsics(1,3)).*inputDepth/cameraIntrinsics(1,1);
camY = (pixY-cameraIntrinsics(2,3)).*inputDepth/cameraIntrinsics(2,2);
camZ = inputDepth;
worldX = cameraPose(1,1) * camX + cameraPose(1,2) * camY + cameraPose(1,3) * camZ + cameraPose(1,4);
worldY = cameraPose(2,1) * camX + cameraPose(2,2) * camY + cameraPose(2,3) * camZ + cameraPose(2,4);
worldZ = cameraPose(3,1) * camX + cameraPose(3,2) * camY + cameraPose(3,3) * camZ + cameraPose(3,4);

end

