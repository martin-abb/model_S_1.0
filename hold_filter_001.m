%
%   hold_filter_001.m
%
%   Martin Krucinski
%
%   2019-04-22  v001    Implementation of "peak-hold" response filter

filt_ord = 1;
Ts      = 1;
a       = 0.8;
A       = [ a ];
B       = [ 1 - a];
C       = [1];
D       = [0];
filt1   = ss(A,B,C,D,Ts);
x0_scalar   = 0.3;
x0_dim      = [ 1 ];

if filt_ord==2,
    temp = filt1 * filt1;
    filt1 = temp;
    x0_dim      = [ 1 ; 1 ];
end

x0      = x0_scalar * x0_dim;

%u = [ zeros(20,1) ; ones(20,1) ; zeros(20,1)  ];
%u = [ zeros(10,1) ; ones(20,1) ; zeros(10,1); ones(10,1) ; zeros(10,1)  ];
%u = [ zeros(5,1) ; ones(10,1) ; zeros(30,1); ones(5,1) ; zeros(40,1)  ];
u = [ 0 0 0 0.5 0.5 0.5 0.5 0.5 0.45 0.45 0.40 0.40 0.35 0.35 0.30 0.30 0.30 ...
    0.35 0.40 0.45 0.50 0.55 1.0 1.0 1.0 1.0 1.0 1.0 0.95 0.97 0.92 0.96 0.99 0.91 0.93 0.96 0.92 0.99 0.95 0.95 ...
    0 0 0 0 0 0 0 0 0 0 0 0 0]';


N = length(u);
t = 0:(N-1);
y = lsim(filt1,u,t, x0);

figure
p1 = stairs(t, [u y]);
set(p1,'LineWidth',3);
title('Normal linear, causal LP filter')

%------------------------------------
%   Reserve space for results

all_y   = zeros(N,N);

%------------------------------------
%   Main filter loop

figure
hold on
p2  = stairs(t, u, 'r');
set(p2,'LineWidth',6);

for i=2:(N-1),
    disp([ 'i = ' num2str(i) ])
    %    i = 35;
    x0  = u(i)*x0_dim;
    upf  = u(i:N);
    tpf  = t(i:N);
    ypf  = lsim(filt1, upf, tpf, x0);
    
    % prepare for reverse filtering
    upr  = flipud(u(1:i));
    tpr  = t(1:i);
    ypr  = lsim(filt1, upr, tpr, x0);
    
    % assemble final simulated output
    
    yprr = flipud(ypr);
    
    y   = [ yprr(1:(i-1)) ; ypf ];
    
    if 0, 
        p3 = stairs(t, y);
    
        set(p3,'LineWidth',3);
    end
    
    %pause
    
    %all_y(i,:)  = y;
    all_y(:,i)  = y;
end

ymax    = max(all_y')';

%figure

%p4 = stairs(t, ymax, 'g');
p4 = plot(t, ymax, 'g');
set(p4,'LineWidth',6);

