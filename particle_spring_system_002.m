%
%   particle_spring_system_002.m
%
%   Martin Krucinski
%
%   2019-04-18  v001    Inital version of discrete-time particle - spring
%   system
%   Start with example of 4 masses
%
%   2019-04-23  v002    2nd version with auto-generation of A-matrix for
%                       systems with arbitrary sizes

Nlip        = 150%100;      % number of points along the lip
Ts          = 0.1;
t_final     = 30;          % simulation end time
Nt          = round(t_final / Ts);           % number of time-points

g           = 9.81;


m           = zeros(Nlip,1);
k           = zeros(Nlip,1);
kv          = zeros(Nlip,1);


for i = 1:Nlip,
    m(i)        = 1;
    k(i)        = 1000;%0;%10;
    kv(i)       = 100;%0.1; %**10
end



%   States are [ ... y_i-1   v_i-1   y_i   v_i   y_i+1   v_i+1 ...]

%   For Nm = 4 (number of masses) this is
%   x = [ y1  v1  y2  v2  y3  v3  y4  v4  ]

N           = Nlip * 2;

x0          = -5/4 * 2 * repmat([1 0], Nlip, 1);

A8(1,:)      = [ 0 1 0 0 0 0 0 0];
A8(3,:)      = [ 0 0 0 1 0 0 0 0];
A8(5,:)      = [ 0 0 0 0 0 1 0 0];
A8(7,:)      = [ 0 0 0 0 0 0 0 1];

A8(2,:)      = [ -k(1)/m(1)  -kv(1)/m(1)     +k(1)/m(1)          0           0                   0           0                   0 ];
A8(4,:)      = [ +k(1)/m(2)  0               (-k(2)-k(1))/m(2)   -kv(2)/m(2) +k(2)/m(2)          0           0                   0 ];
A8(6,:)      = [ 0           0               +k(2)/m(3)          0           (-k(3)-k(2))/m(3)   -kv(3)/m(3) +k(3)/m(3)          0 ];
%A8(8,:)      = [ 0           0               0                   0           +k(3)/m(4)          0           (-k(4)-k(3))/m(4)   -kv(4)/m(4)];
A8(8,:)      = [ 0           0               0                   0           +k(3)/m(4)          0           (-k(3))/m(4)   -kv(4)/m(4)];

A               = zeros(N,N);

for i = 1:Nlip,
    
    ix      = 1 + (i-1)*2;
%    ic      = (ix-1);
    %   odd ix
    A(ix,ix:(ix+1))  = [  0 1 ];
    A(ix+1,ix+1)    = -kv(i) / m(i);
    %   even ix
    if ix==1,
        A(ix+1,1)   = -k(1)/m(1);
        A(ix+1,3)   = +k(1)/m(1);
    elseif ix==(N-1),
        A(ix+1, ix-2)   = +k(i-1)/m(i);
        A(ix+1, ix)     = (-k(i-1))/m(i);
    else % ix == >1 && ix < (N-1)
        A(ix+1, ix-2)   = +k(i-1)/m(i);
        A(ix+1, ix)     = (-k(i)-k(i-1))/m(i);
        A(ix+1, ix+2)   = +k(i)/m(i);
        
    end
    
end



%   Input Matrix, driving acceleration is -g

%B           = [ 0 1 0 1 0 1 0 1]' * -g;
B          = -g * repmat([0 ; 1], Nlip, 1);

%   Output matrix

%C           = eye(N);
% C           = [     1 0 0 0  0 0 0 0 ; ...
%     0 0 1 0  0 0 0 0 ; ...
%     0 0 0 0  1 0 0 0 ; ...
%     0 0 0 0  0 0 1 0 ];

C   = zeros(Nlip, N);

for i = 1:Nlip,
    
    ix      = 1 + (i-1)*2;
    C(i,ix)     = 1;
end


%   Direct feed-through matrix

%D           = zeros(N,1);
D           = zeros(Nlip,1);

sysC        = ss(A, B, C, D);
sysC.InputName  = 'g Enable';

%sysC.StateName = { 'y1', 'v1', 'y2', 'v2', 'y3', 'v3', 'y4', 'v4' };

%sysC.OutputName = sysC.StateName;
%sysC.OutputName = { 'y1', 'y2', 'y3', 'y4' };


sysD        = c2d(sysC, Ts);

%   Define object the lip will land on


%obj       = [2 1 1.4 1.3]';
%obj       = [10 0 0 10]';
obj         = [1 ; zeros(Nlip-2,1) ; 1];
obj(15) = 0.8
obj(32) = 0.95

objmax    = max(obj);
lip0    = objmax * 1.02; %***1.1;
%x0      = [lip0 0 lip0 0 lip0 0 lip0 0]';
x0      = lip0 * repmat([1 ; 0], Nlip, 1);

x       = x0;
figure
hold on
plot(obj, 'r');
all_y   = [];

all_t = Ts*(0:(Nt-1))';

disp('Starting simulation...')

for t=1:Nt,
    disp(all_t(t))
    y       = sysD.C * x;
    plot(y)
%    xprev   = x;
    xnext   = sysD.A * x + sysD.B * 1;
    ynext   = sysD.C * xnext;
    for i=1:Nlip,
        ix      = 1 + (i-1)*2;
        if xnext(ix) < obj(i),
            xnext(ix) = obj(i);         % constraint mass particle position to top of object
            xnext(ix+1)     = 0;        % set mass particle velocity to zero
        end
    end

    x = xnext;
    
    all_y(:,t) = y;
    
 %   pause
end

    
figure
stairs(all_t, all_y')
xlabel('t')
ylabel('y [m]')
title('Suction cup lip mass particle positions [m]')


    