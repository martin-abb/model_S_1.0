function [PickPoints_Bin,rows]  = Filter_Bin_Pick_Points(PickPoints, Bin)
%function [PickPoints_Bin,pick_set_X_Bin,pick_set_Y_Bin]  = Filter_Bin_Pick_Points(PickPoints, Bin)
%Filter_Bin_Pick_Points Remove all pick points that are outside the bin
%(within tolerances Bin.dx, dy, dz)

X_min                 = Bin.X - Bin.dx;
X_max                 = Bin.X + Bin.BW + Bin.dx;

Y_min                 = Bin.Y - Bin.dy;
Y_max                 = Bin.Y + Bin.BL + Bin.dy;

Z_min                 = Bin.Z - Bin.dz;
Z_max                 = Bin.Z + Bin.BH + Bin.dz;


PX                  = PickPoints(:,1);
PY                  = PickPoints(:,2);
PZ                  = PickPoints(:,3);

%     [rows, cols]        = find( (PX >= X_min) & (PX <= X_max) & ...
%                                 (PY >= Y_min) & (PY <= Y_max) & ...
%                                 (PZ >= Z_min) & (PZ <= Z_max) );
%
%
%     % The syntax below DOES NOT WORK, RETURNS SQUARE matrices with N*N
%     % elements..
%     % lipX                = pickX(rows, cols);
%     % lipY                = pickY(rows, cols);
%     % lipZ                = pickZ(rows, cols);
%
%     %   instead try
%
%     S           = size(PickPoints);
%     linearInd   = sub2ind(S, rows, cols);
%     FPX        = PX(linearInd);
%     FPY        = PY(linearInd);
%     FPZ        = PZ(linearInd);
%

rows = find( ...
    (PX >= X_min) & (PX <= X_max) & ...
    (PY >= Y_min) & (PY <= Y_max) & ...
    (PZ >= Z_min) & (PZ <= Z_max) );


% The syntax below DOES NOT WORK, RETURNS SQUARE matrices with N*N
% elements..
% lipX                = pickX(rows, cols);
% lipY                = pickY(rows, cols);
% lipZ                = pickZ(rows, cols);

%   instead try

% S           = size(PickPoints);
% linearInd   = sub2ind(S, rows, cols);
FPX        = PX(rows);
FPY        = PY(rows);
FPZ        = PZ(rows);

PickPoints_Bin  = [ FPX FPY FPZ ];

end

