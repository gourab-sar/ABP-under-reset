clear all;
clc;
format long
clear variables;

tic

h = 0.01;           
  
% resetting rate
R = 2;

% threshold error
eps = 0.1;


v = 1.0; % velocity
DR = 4; %rotational diffusion constant
a = 5.0; % length of physical space

tt = [];

m = 0;
x = [];
y = [];
theta = [];
j = 1;

% initial conditions
x(j) = 0;
y(j) = 0;
theta(j) = -pi + 2*pi*rand();
%theta(j) = 0;

while 1
    T = round(exprnd(R),2)/h;
    if (T>=1)
        break
    end        
end
% T2 = [T];
% time iteration
while j<=100000
    x(j+1) = x(j) + v*cos(theta(j))*h;
    y(j+1) = y(j) + v*sin(theta(j))*h;
    theta(j+1) = theta(j) + sqrt(2*DR*h)*randn(1,1);

    m = m+1;
    
    if (m == round(T))
        x(j+1) = 0;
        y(j+1) = 0;
        theta(j+1) = -pi + 2*pi*rand();
        %theta(j+1) = 0;
        while 1
            T1 = round(exprnd(R),2)/h;
            if (T1>=1)
                break
            end        
        end
%         T2 = [T2 T1];
        T = T + T1;
    end
    
    % Boundary conditions: box vertical edge
    if (abs(x(j+1))>a)
        x(j+1) = sign(x(j+1))*(2*a-abs(x(j+1)));
    end
    % Boundary conditions: box horizontal edge
    if (abs(y(j+1))>a)
        y(j+1) = sign(y(j+1))*(2*a-abs(y(j+1)));
    end
   
    j = j+1;
    
end
time = m*h; %first passage time

data = [[0:1:j-1]'.*h x' y'];
nData = size(data,1); %# number of data points
numberOfdeltaT = floor(nData/4); %# for MSD, dt should be up to 1/4 of number of data points

msd = zeros(numberOfdeltaT,5); %# We'll store [mean, std, n]

%# calculate msd for all deltaT's

for dt = 1:numberOfdeltaT
   deltaCoords = data(1+dt:end,2:3) - data(1:end-dt,2:3);
   squaredDisplacement = deltaCoords.^2;
   sumsquaredDisplacement = sum(squaredDisplacement,2); %# dx^2+dy^2+dz^2
   

   msd(dt,1) = mean(sumsquaredDisplacement); %# average
   msd(dt,2) = std(sumsquaredDisplacement); %# std
   msd(dt,3) = length(sumsquaredDisplacement); %# n
   msd(dt,4) = mean(squaredDisplacement(:,1)); %x msd
   msd(dt,5) = mean(squaredDisplacement(:,2)); %y msd
end

% figure(1)
% % dlmwrite('fpt_expreset_pt1_BC1_protocol1.dat',tt,'delimiter','\t')
% histogram(sqrt(x.^2+y.^2),'Normalization', 'probability')
% xlabel('d')
% ylabel('\it{p(r)} [n.u.]')
% title(['Histogram plot for R = ' num2str(R)])

% figure(2)
% dlmwrite('fpt_expreset_pt1_BC1_protocol1.dat',tt,'delimiter','\t')
% histogram(sqrt(x.^2+y.^2),'Normalization', 'pdf')
% xlabel('d')
% ylabel('\it{p(r)} [n.u.]')
% title(['Histogram plot for R = ' num2str(R)])
% figure
% plot(x,y,'linewidth',2)
% hold on
% plot(0,0,'o')
% axis([-a-0.01 a+0.01 -a-0.01 a+0.01]);
toc