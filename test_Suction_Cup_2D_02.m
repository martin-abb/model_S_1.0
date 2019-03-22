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

%--------------------------------------------------------------------------------------
%   Plot object surface, top view
f1          = figure;
image(object1*64)
colormap gray


lip1 = imgaussfilt(object1, N_Gauss);

%--------------------------------------------------------------------------------------
%   plot 1D plot of object profile and suction cup lip profile
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

%--------------------------------------------------------------------------------------
%   Calculate & overlay magnitude of leakage distance
d1 = abs(object1 - lip1);
figure(f2)
plot(d1,'g')

%--------------------------------------------------------------------------------------
%   plot leakage magnitude surface, top view
f4          = figure;
image(d1*64)
colormap gray
ax=gca;
set(ax, 'Ydir', 'normal');

%--------------------------------------------------------------------------------------
%   Define Object 2

temp1   = [ the_ones the_zeros the_ones the_zeros the_ones the_zeros the_ones ];
temp2   = [ the_ones the_ones  the_ones the_zeros the_ones the_zeros the_ones ];
temp3   = [ the_ones the_ones  the_ones the_ones  the_ones the_zeros the_ones ];

temp4   = ( temp1' * the_ones )';
temp5   = ( temp2' * the_ones )';
temp6   = ( temp3' * the_ones )';

object2 = [ temp4 ; temp5 ; temp6 ];

%--------------------------------------------------------------------------------------
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

%--------------------------------------------------------------------------------------
%   Define Object 4, circular donut with an internal leakage path, located
%   between 0 and 90 degrees

deg         = pi/180;
mm          = 0.001;

init_Suction
N4          = N*20;     % set size of spatial grid
object4     = zeros(N4,N4);
Lx          = Suction.RO * 2 * 1.5; %   define linear dimension of modeled region
%   to be 150% of suction cup outer
%   diameter
Ly          = Lx;


%   fill in area of suction cup lip with height 1
%   and fill in leakage channel

plug_channel    = true;     % true - plug up channel so NO leakage occurs, false - leave channel open

R1      = Suction.RI + 1/3*(Suction.RO - Suction.RI);   % inner radius of leakage channel
R2      = Suction.RI + 2/3*(Suction.RO - Suction.RI);   % outer radius of leakage channel
th_leak = 5*deg;        % angular size of leakage channels at top and right side

for x=1:N4,
    for y=1:N4,
        xl  = (x/N4-0.5)*Lx;      % relative linear x -coordinate
        yl  = (y/N4-0.5)*Ly;
        Rl   = sqrt(xl^2+yl^2);
        if ( (Rl >= Suction.RI) && (Rl <= Suction.RO) ),
            object4(y,x)    = 1;
        end
        
        %   add leakage channel
        
        theta   = atan2(yl,xl);    % (x,y) angle from 0 to 2*pi
        if theta<0,
            theta   = 2*pi + theta; % convert the 0 to -180 range to 360 to 180
        end
        if ( (Rl >= R1) && (Rl <= Suction.RO) && (theta > (90*deg - th_leak) && (theta < (90*deg + th_leak) ) ) ),   % upper leakage path
            object4(y,x)    = 0;
        elseif ( (Rl >= R1) && (Rl <= R2) && (theta > 0*deg) && (theta < 90*deg) ),   % radial leakage path
            object4(y,x)    = 0;
        elseif ( (Rl >= Suction.RI) && (Rl <= R2) && ( (theta > ( 360*deg - th_leak) || (theta < th_leak) ) ) ),   % right side leakage path
            object4(y,x)    = 0;
            
        end
        
        if plug_channel,
            if ( (Rl >= R1) && (Rl <= R2) && (theta > 40*deg) && (theta < 50*deg) ),   % radial leakage path
                object4(y,x)    = 1;
            end
        end
        
    end
end




%   Override Object 2 for tests
disp('Overriding object2 with object4 !')
object2     = object4;

%disp('Overriding object2 with RANDOM object !')
%object2     = rand(size(object3));


lip2 = imgaussfilt(object2, N_Gauss);
d2 = abs(object2 - lip2);


%--------------------------------------------------------------------------------------
%   plot top view of object
f5          = figure;
image(object2*64)
colormap gray
axis equal
ax=gca;
set(ax, 'Ydir', 'normal');

%--------------------------------------------------------------------------------------
%   plot top view of suction cup surface
f6          = figure;
image(lip2*64)
colormap gray
axis equal
ax=gca;
set(ax, 'Ydir', 'normal');

%--------------------------------------------------------------------------------------
%   plot top view of suction cup leakage path height amplitude
f7          = figure;
image(d2*64)
colormap gray
axis equal
ax=gca;
set(ax, 'Ydir', 'normal');

%   generate BINARY difference image

leak_threshold  = 0.2%0.3;
bd2         = d2 > leak_threshold;

%--------------------------------------------------------------------------------------
%   plot top view of THRESHOLDED leakage path height
f8          = figure;
image(bd2*64/2)
colormap gray
axis equal
ax=gca;
set(ax, 'Ydir', 'normal');

%--------------------------------------------------------------------------------------
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

connectedness = 2;  % 1 - rectangular , 2 - circular suction cup

switch connectedness,
    
    case 1,
        %--------------------------------------------------------------------------------------
        %   RECTANGULAR CONNECTEDNESS
        %   Find regions that connect OUTSIDE (i.e. TOP of image)
        %   with INSIDE (i.e. BOTTOM of image)
        
        CC2     = bwconncomp(bd2,4);
        CC2
        
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
        
    case 2,
        %--------------------------------------------------------------------------------------
        %   CIRCULAR SUCTION CUP  CONNECTEDNESS
        %   Find regions that connect OUTSIDE (i.e. TOP of image)
        %   with INSIDE (i.e. BOTTOM of image)
        
        %   threshold the area outside the suction cup to white (i.e. LEAKAGE area)
        
        for x=1:N4,
            for y=1:N4,
                xl  = (x/N4-0.5)*Lx;      % relative linear x -coordinate
                yl  = (y/N4-0.5)*Ly;
                Rl   = sqrt(xl^2+yl^2);
                if ( (Rl >= Suction.RO) || (Rl <= Suction.RI) ),
                    bd2(y,x)    = true;
                end
            end
        end
        
        
        
        CC2     = bwconncomp(bd2,4);
        CC2
        
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
            
            iXX = find(N4/2==XX);
            iYY = find(N4/2==YY);
            intersect_iXX_iYY = intersect(iXX, iYY);
            
            if ( (min(YY)==1) && (length(intersect_iXX_iYY)>0) ), % Connected path from OUTSIDE to INSIDE exists
                disp('Found connecting leakage path!');
                num_leak        = num_leak + 1;
                all_leak_no(num_leak)   = r;
                all_leak_Idx{num_leak}  = Idx_region;
                all_leak_X{num_leak}    = XX;
                all_leak_Y{num_leak}    = YY;
            end
        end
        
        disp([ 'Found ' num2str(num_leak) ' leakage paths' ] );
        
end % switch connectedness

%--------------------------------------------------------------------------------------
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
axis equal
ax=gca;
set(ax, 'Ydir', 'normal');




