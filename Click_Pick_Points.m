function [PickPoints,pick_set_X,pick_set_Y]    = Click_Pick_Points(worldX, worldY, worldZ)
%Click_Pick_Points Select pick points by clicking in the RGB image
%   Detailed explanation goes here

disp('Click on desired pick points, press ENTER when done...')

gx          = ginput;
pick_set_X  = floor(gx(:,1));
pick_set_Y  = floor(gx(:,2));

hold on
set(gcf,'DefaultLineLineWidth',3);
plot(pick_set_X, pick_set_Y, 'ro')

S           = size(worldX);

N_pick      = length(pick_set_X);

linearInd   = sub2ind(S, pick_set_Y, pick_set_X);

pickX       = reshape(worldX(linearInd),N_pick,1);
pickY       = reshape(worldY(linearInd),N_pick,1);
pickZ       = reshape(worldZ(linearInd),N_pick,1);

PickPoints  = [ pickX pickY pickZ];

end

