%
%   test8.m
%
%   Testing of Suction 1.0
%   Martin Krucinski
%
%   2019-01-23      v8

img_raw_color   = imread('24_image_raw.color.png');
img_raw_depth   = imread('24_image_raw.depth.png');
img_proc_pick   = imread('24_image_proc.pickpoint.png');
scale           = 1/10^4;
inputDepth      = double(img_raw_depth)*scale;

cameraIntrinsics    = load('test-camera-intrinsics.txt','ascii');
cameraPose          = load('test-camera-pose.txt','ascii');

N_metrics           = 4;    % number of metrics to record for each pick point

Suction_Cup_Check = true;

if Suction_Cup_Check,
    disp('INCLUDING check for suction cup lip completeness...')
else
    disp('SKIPPING check for suction cup lip completeness...')
end

%----------------------------------------------------------------------
f1=figure;
image(img_raw_color);
title('Raw Color Image')

%----------------------------------------------------------------------
f2=figure;
image(img_raw_depth/256);
colormap('gray')
title('Raw Depth Image')

%----------------------------------------------------------------------
%
%   Filter depth image to move glitch points at d = 0.0 further back, to
%   e.g. 2.0 m
%
%   Filter the inputDepth image, not the img_raw_depth


d_move_back     = 20%8.0;   %   Move back distance
inputDepth_filtered     = filter_depth( inputDepth, d_move_back);

%img_depth_filtered = inputDepth_filtered / scale;      % create filtered depth image map
img_filtered_depth = inputDepth_filtered / scale;      % create filtered depth image map

%----------------------------------------------------------------------
f3=figure;
image(img_filtered_depth/256);
colormap('gray')
title('Filtered Depth Image (glitch filtered)')

%----------------------------------------------------------------------
f4=figure;
image(img_proc_pick);
title('MIT Princeton Pick Points')

