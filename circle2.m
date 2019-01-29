function [Cx,Cy] = circle2(R, dtheta)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here
   theta               = [0:dtheta:2*pi];
    Cx                 = R*cos(theta);
    Cy                 = R*sin(theta);
end

