clear all;
clc;
format long
clear variables;

tic

h=0.01;
  
% target
% xt = 4.0;
% yt = 4.0;

%threshold error
% eps = 0.1;


v = 1.0; % velocity
DR = 1; %rotational diffusion constant
a = 5.0; % length of physical space

% initial conditions
x(1)=0;
y(1)=0;
theta(1)=-pi+2*pi*rand();
j = 1;
m = 0;
% time iteration
while j <= 50000
    x(j+1)=x(j)+v*cos(theta(j))*h;
    y(j+1)=y(j)+v*sin(theta(j))*h;
    theta(j+1)=theta(j)+sqrt(2*DR*h)*randn(1,1);%wgn(1,1,0);
    m = m+1; 
    
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
% plot(x,y,'linewidth',2)

% XX = [];
% 
% for i = 1:j
%     if (abs(y(i)) < 0.1)
%         XX = [XX x(i)];
%     end
% end
data = [[0:1:j-1]'.*h x' y'];
% figure(1)
% histogram(sqrt(x.^2+y.^2),'Normalization', 'probability')
% xlabel('x')
% ylabel('\it{p(r)} [n.u.]')

% figure(2)
% histogram(sqrt(x.^2+y.^2),'Normalization', 'pdf')
% xlabel('x')
% ylabel('\it{p(r)} [n.u.]')
% title(['Histogram plot for V = ' num2str(v)])
toc

nData = size(data,1); %# number of data points
numberOfdeltaT = floor(nData/4); %# for MSD, dt should be up to 1/4 of number of data points

msd = zeros(numberOfdeltaT,3); %# We'll store [mean, std, n]

%# calculate msd for all deltaT's

for dt = 1:numberOfdeltaT
   deltaCoords = data(1+dt:end,2:3) - data(1:end-dt,2:3);
   squaredDisplacement = sum(deltaCoords.^2,2); %# dx^2+dy^2+dz^2

   msd(dt,1) = mean(squaredDisplacement); %# average
   msd(dt,2) = std(squaredDisplacement); %# std
   msd(dt,3) = length(squaredDisplacement); %# n
end