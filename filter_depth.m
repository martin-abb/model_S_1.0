function [filteredDepth] = filter_depth(inputDepth,d_move_back)
%filter_depth Filter depth image to move glitch points at d = 0.0 further back, to d_move_back

izeros = find(inputDepth == 0);      % find glitches in depth data

filteredDepth           = inputDepth;
filteredDepth(izeros)   = d_move_back*ones(size(izeros));
filteredDepth           = reshape(filteredDepth, size(inputDepth));

%***disp('FILTERING Depth Map!!!')

