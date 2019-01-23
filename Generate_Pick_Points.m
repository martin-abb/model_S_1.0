function [pickPoints,pick_set_X,pick_set_Y] = Generate_Pick_Points(NskipX, NskipY, worldX, worldY, worldZ)
% Generate_Pick_Points for heat maps every NskipX and NskpY pixels in the
% raw image
%   worldX etc. assumed to be same dimensions as original RGB and D images

[NY,NX]=size(worldX);

pick_set_X = 1:NskipX:NX;      % indicies in image x axis to use for pick point grid
pick_set_Y = 1:NskipY:NY;

N_pick      = length(pick_set_X) * length(pick_set_Y);

pickX       = reshape(worldX(pick_set_Y, pick_set_X),N_pick,1);
pickY       = reshape(worldY(pick_set_Y, pick_set_X),N_pick,1);
pickZ       = reshape(worldZ(pick_set_Y, pick_set_X),N_pick,1);

pickPoints  = [ pickX pickY pickZ];


end

