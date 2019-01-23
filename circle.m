function [Cx,Cy] = circle(R)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here
   theta               = [0:0.01:2*pi];
    Cx                 = R*cos(theta);
    Cy                 = R*sin(theta);
end

