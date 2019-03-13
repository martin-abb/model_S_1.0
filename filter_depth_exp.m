function [filteredDepth] = filter_depth_exp(inputDepth,d_move_back)
%filter_depth_exp Filter depth image in MY experimental data
%to move glitch points at d = x_min = 0.2 further back, to d_move_back

izeros = find(inputDepth < 0.81);      % find glitches in depth data

filteredDepth           = inputDepth;
filteredDepth(izeros)   = d_move_back*ones(size(izeros));
filteredDepth           = reshape(filteredDepth, size(inputDepth));

%***disp('FILTERING Depth Map!!!')

