clear all
x = load('fpt_reset_data.dat');
R = x(:,1);
x1 = load('fpt_data.dat');
fpt = zeros(length(R),1);
h = 0.01;
for i = 1:length(R)
    y3 = exprnd(R(i),2*length(x1),1);
    y2 = y3(y3>=h);
    y1 = y2(1:length(x1),1);
%     y1 = exprnd(R(i),length(x1),1);
    A = min(x1,y1);
    B = x1<y1;
    fpt(i) = mean(A) / (sum(B)/length(x1));
end
figure()
loglog(x(:,1),x(:,2),'-o','Linewidth',2,'Markersize',8)
hold on
loglog(x(:,1),fpt,'-*','Linewidth',2,'Markersize',8)
set(gca,'Linewidth',2,'Fontsize',18)
