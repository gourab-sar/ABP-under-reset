clc;
format long
clear variables;

tic

h = 0.01;  

r = 0.25;
t2 = [];
for l=1:1
ttt = 0;
l
for k=1:10
    k
    angle(k) = 2*pi*k/10;
% target
xt = r*cos(angle(k));
yt = r*sin(angle(k));

% threshold error
eps = 0.1;


v = 1.0; % velocity
DR = 1.0; %rotational diffusion constant
a = 5.0; % length of physical space

% number of realisations
realisation = 200;

time = zeros(realisation,1);
tt = [];


for i = 1:realisation
i;

m = 0;
x = [];
y = [];
theta = [];
j = 1;

% initial conditions
x(j) = 0.0;
y(j) = 0.0;
theta(j) = -pi + 2*pi*rand();
% time iteration
while 1   
    x(j+1) = x(j) + v*cos(theta(j))*h;
    y(j+1) = y(j) + v*sin(theta(j))*h;
    theta(j+1) = theta(j) + sqrt(2*DR*h)*randn(1,1);

    m = m+1;
    
    % Boundary conditions: box vertical edge
    if (abs(x(j+1))>a)
        x(j+1) = sign(x(j+1))*(2*a-abs(x(j+1)));
    end
    % Boundary conditions: box horizontal edge
    if (abs(y(j+1))>a)
        y(j+1) = sign(y(j+1))*(2*a-abs(y(j+1)));
    end
   
    distance = norm([xt,yt]-[x(j+1),y(j+1)]); % distance from the target
    
    if (distance < eps)
        break;
    end
    
    j = j+1;
    
end
time(i) = m*h; %first passage time
end
tt = [mean(time) std(time) std(time)/mean(time)];

ttt = ttt+tt;
end
ttt = ttt/10;
t2 = [t2;r,ttt];
r = r+0.25;
end
dlmwrite('dist_vs_fpt_bc2.dat',t2,'delimiter','\t')
toc