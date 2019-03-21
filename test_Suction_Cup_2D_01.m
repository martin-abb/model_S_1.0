%
%   test_Suction_Cup_2D_01.m
%
%   Martin Krucinski
%
%

N=20;
N_Gauss     = N/4;

the_zeros   = zeros(1,N);
the_ones    = ones(1,N);
object1  = [the_zeros the_ones the_zeros the_ones the_zeros];
f1          = figure;
image(object1*64)
colormap gray


lip1 = imgaussfilt(object1, N_Gauss);

f2          = figure;
set(gcf, 'DefaultLineLineWidth',3);
plot(object1,'b')
hold on
plot(lip1,'r')

% object1_partial = object1(1:(N+1));
% lip1_partial = imgaussfilt(object1_partial, 1);
% 
% f3          = figure;
% set(gcf, 'DefaultLineLineWidth',3);
% hold on
% plot(object1_partial,'b')
% plot(lip1_partial,'r')

d1 = abs(object1 - lip1);
figure(f2)
plot(d1,'g')

f4          = figure;
image(d1*64)
colormap gray

%   Define Object 2

temp1   = [ the_ones the_zeros the_ones the_zeros the_ones the_zeros the_ones ];
temp2   = [ the_ones the_ones  the_ones the_zeros the_ones the_zeros the_ones ]; 
temp3   = [ the_ones the_ones  the_ones the_ones  the_ones the_zeros the_ones ]; 

temp4   = ( temp1' * the_ones )';
temp5   = ( temp2' * the_ones )';
temp6   = ( temp3' * the_ones )';

object2 = [ temp4 ; temp5 ; temp6 ];

%   Define Object 3

t0      = the_zeros;
t1      = the_ones;
temp1   = [ t1 t0 t1 t1 t1 t1 t1 t1 t1 t0 t1];
temp2   = [ t1 t0 t0 t1 t1 t0 t0 t0 t1 t0 t1];
temp3   = [ t1 t1 t0 t0 t0 t0 t1 t0 t1 t0 t1];
temp4   = [ t1 t1 t1 t1 t1 t1 t1 t0 t1 t1 t1];

temp5   = ( temp1' * t1 )';
temp6   = ( temp2' * t1 )';
temp7   = ( temp3' * t1 )';
temp8   = ( temp4' * t1 )';


object3 = [ temp5 ; temp6 ; temp7 ; temp8 ];

%   Override Object 2 for tests
% disp('Overriding object2 with object3 !')
% object2     = object3;

%disp('Overriding object2 with RANDOM object !')
%object2     = rand(size(object3));


lip2 = imgaussfilt(object2, N_Gauss);
d2 = abs(object2 - lip2);


f5          = figure;
image(object2*64)
colormap gray


f6          = figure;
image(lip2*64)
colormap gray


f7          = figure;
image(d2*64)
colormap gray

%   generate BINARY difference image

leak_threshold  = 0.2%0.3;
bd2         = d2 > leak_threshold;
f8          = figure;
image(bd2*64/2)
colormap gray

%   Find connected regions

SZ                  = size(object2);

%rangeX              = 0:dx:(L);
%rangeY              = 0:dx:(W);

rangeX          = 1:SZ(2);
rangeY          = 1:SZ(1);

% plot_X_range = [ rangeX(1) - 0.1 , rangeX(end) + 0.1 ];
% plot_Y_range = [ rangeY(1) - 0.1 , rangeY(end) + 0.1 ];

NX                  = length(rangeX);
NY                  = length(rangeY);

[WX,WY]     = meshgrid( rangeX , rangeY );

CC2     = bwconncomp(bd2,4);
CC2

%   Find regions that connect OUTSIDE (i.e. TOP of image)
%   with INSIDE (i.e. BOTTOM of image)

Nregions    = CC2.NumObjects;

all_leak_Idx = {};
all_leak_X  = {};
all_leak_Y  = {};
num_leak    = 0;
all_leak_no = [];

for r = 1:Nregions,
    Idx_region  = CC2.PixelIdxList{r};
    XX  = WX(Idx_region);
    YY  = WY(Idx_region);
    
    if ( (min(YY)==1) && (max(YY)==NY) ), % Connected path from TOP to BOTTOM exists
        disp('Found connecting leakage path!');
        num_leak        = num_leak + 1;
        all_leak_no(num_leak)   = r;
        all_leak_Idx{num_leak}  = Idx_region;
        all_leak_X{num_leak}    = XX;
        all_leak_Y{num_leak}    = YY;
    end
end

disp([ 'Found ' num2str(num_leak) ' leakage paths' ] );

f9 = figure;

%***imgleaks = zeros(size(bd2));
imgleaks = double(bd2)*0.5;


for l   = 1:num_leak,
    leak_Idx        = all_leak_Idx{l};
    SLI             = size(leak_Idx);
    imgleaks(leak_Idx)  = ones(SLI) * 1.0;
end

image(imgleaks*64)
colormap gray



        
