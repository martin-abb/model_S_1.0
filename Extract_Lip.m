function [ lipX, lipY, lipZ ] = Extract_Lip( Pick_PC_Rel_X, Pick_PC_Rel_Y, Pick_PC_Rel_Z, Suction)
%Extract_Lip Extract the point cloud under the lip of the suction cup

R                   = sqrt(Pick_PC_Rel_X.^2 + Pick_PC_Rel_Y.^2);
[rows, cols]        = find( (R>=Suction.RI) & (R<=Suction.RO) );

S           = size(Pick_PC_Rel_X);
linearInd   = sub2ind(S, rows, cols);
lipX        = Pick_PC_Rel_X(linearInd);
lipY        = Pick_PC_Rel_Y(linearInd);
lipZ        = Pick_PC_Rel_Z(linearInd);

end

