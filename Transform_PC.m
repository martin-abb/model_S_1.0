function [Pick_PC_Rel_X, Pick_PC_Rel_Y, Pick_PC_Rel_Z ] = Transform_PC(WorldX, WorldY, WorldZ, PickPose)
%Transform_PC Transform input point cloud to a relative point cloud around
%the pick point and approach angle specified in PickPose
%   Detailed explanation goes here
Pick_PC_Rel_X = PickPose(1,1) * worldX + PickPose(1,2) * worldY + PickPose(1,3) * worldZ + PickPose(1,4);
Pick_PC_Rel_Y = PickPose(2,1) * worldX + PickPose(2,2) * worldY + PickPose(2,3) * worldZ + PickPose(2,4);
Pick_PC_Rel_Z = PickPose(3,1) * worldX + PickPose(3,2) * worldY + PickPose(3,3) * worldZ + PickPose(3,4);
end

