%
%   particle_spring_system_001.m
%
%   Martin Krucinski
%
%   2019-04-18  v001    Inital version of discrete-time particle - spring
%   system

%   Start with example of 4 masses

m(1)        = 1;
k(1)        = 10;
kv(1)       = 1;

m(2)        = 1;
k(2)        = 10;
kv(2)       = 1;

m(3)        = 1;
k(3)        = 10;
kv(3)       = 1;

m(4)        = 1;
k(4)        = 10;
kv(4)       = 1;

%kv = [0 0 0 0]
%k  = [0 0 0 0]

g           = 9.81;
Ts          = 0.01;



%   States are [ ... y_i-1   v_i-1   y_i   v_i   y_i+1   v_i+1 ...]

%   For Nm = 4 (number of masses) this is
%   x = [ y1  v1  y2  v2  y3  v3  y4  v4  ]

x0          = -5/4*[ 2 0 2 0 2 0 2 0]';

A(1,:)      = [ 0 1 0 0 0 0 0 0];
A(3,:)      = [ 0 0 0 1 0 0 0 0];
A(5,:)      = [ 0 0 0 0 0 1 0 0];
A(7,:)      = [ 0 0 0 0 0 0 0 1];

A(2,:)      = [ -k(1)/m(1)  -kv(1)/m(1)     +k(1)/m(1)          0           0                   0           0                   0 ];
A(4,:)      = [ +k(1)/m(2)  0               (-k(2)-k(1))/m(2)   -kv(2)/m(2) +k(2)/m(2)          0           0                   0 ];
A(6,:)      = [ 0           0               +k(2)/m(3)          0           (-k(3)-k(2))/m(3)   -kv(3)/m(3) +k(3)/m(3)          0 ];
%A(8,:)      = [ 0           0               0                   0           +k(3)/m(4)          0           (-k(4)-k(3))/m(4)   -kv(4)/m(4)];
A(8,:)      = [ 0           0               0                   0           +k(3)/m(4)          0           (-k(3))/m(4)   -kv(4)/m(4)];

N           = length(A);

%   Input Matrix, driving acceleration is -g

B           = [ 0 1 0 1 0 1 0 1]' * -g;

%   Output matrix

%C           = eye(N);
C           = [     1 0 0 0  0 0 0 0 ; ...
    0 0 1 0  0 0 0 0 ; ...
    0 0 0 0  1 0 0 0 ; ...
    0 0 0 0  0 0 1 0 ];

%   Direct feed-through matrix

%D           = zeros(N,1);
D           = zeros(4,1);

sysC        = ss(A, B, C, D);
sysC.InputName  = 'g Enable';

sysC.StateName = { 'y1', 'v1', 'y2', 'v2', 'y3', 'v3', 'y4', 'v4' };

%sysC.OutputName = sysC.StateName;
sysC.OutputName = { 'y1', 'y2', 'y3', 'y4' };


sysD        = c2d(sysC, Ts);

%   Define lip


obj       = [2 1 1.4 1.3]';

objmax    = max(obj);
lip0    = objmax*1.1;
x0      = [lip0 0 lip0 0 lip0 0 lip0 0]';
x       = x0;
figure
hold on
plot(obj, 'r');

for t=1:50000,
    disp(t)
    y       = sysD.C * x;
    plot(y)
%    xprev   = x;
    xnext   = sysD.A * x + sysD.B * 1
    ynext   = sysD.C * xnext;
    for i=1:4,
        ix      = 1 + (i-1)*2;
        if obj(i) > xnext(ix),
            xnext(ix) = obj(i);         % constraint mass particle position to top of object
            xnext(ix+1)     = 0;        % set mass particle velocity to zero
        end
    end

    x = xnext;
    pause
end

    
    