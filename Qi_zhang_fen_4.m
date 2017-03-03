% function [T, Y]=deng_fen_szJ(parameters, orders, TSim, Y0)
%% This program is to solve numerical solutions of Fractional differential equation for paper of LingMing Deng,
%% This program  is completed on Apr. 27, 2011, by Dong Li
 
function [T, Y]=Qi_zhang_fen_4( )
parameters=[14 43 -1 16 4]; %parameters of the system:
orders=[0.95 0.95 0.95];%orders of derivatives
% orders=[1 1 1];
TSim=30;
Y0=[0.1  0.2  0.3];
%
% Numerical Solution of the Fractional differential System
%
% A=[ -a      a         0;
%     c       d         0;
%     0      0         -b];
% 
% Fai = [r*X(2)*X(3)  -X(1)*X(3)  X(1)*X(2)]';
% y = A*X + Fai;
%

% time step:
h=0.001; 
% number of calculated mesh points:
n=round(TSim/h);
%orders of derivatives, respectively:
q1=orders(1); 
q2=orders(2); 
q3=orders(3);
% constants of the system:
a=parameters(1); 
b=parameters(2); 
c=parameters(3);
d=parameters(4);
r=parameters(5);
 
% binomial coefficients calculation:
cp1=1; cp2=1; cp3=1;
for j=1:n
    c1(j)=(1-(1+q1)/j)*cp1;
    c2(j)=(1-(1+q2)/j)*cp2;
    c3(j)=(1-(1+q3)/j)*cp3;
  
    cp1=c1(j); cp2=c2(j); cp3=c3(j);
end
% initial conditions setting:
x(1)=Y0(1); 
y(1)=Y0(2); 
z(1)=Y0(3);
% calculation of phase portraits /numerical solution/:
for i=2:n

    x(i)=(-a*x(i-1)+a*y(i-1)+r*y(i-1)*z(i-1))*h^q1 - memo(x, c1, i);
    y(i)=(c*x(i)+d*y(i-1)-x(i)*z(i-1))*h^q2 - memo(y, c2, i);
    z(i)=(-b*z(i-1)+x(i)*y(i))*h^q3 - memo(z, c3, i);
    
end
for j=1:n
    Y(j,1)=x(j);
    Y(j,2)=y(j);
    Y(j,3)=z(j);
end
% T=h:h:TSim;
subplot(3,1,1);
plot(Y(:,1),'b');
xlabel('t');
ylabel('x_{1}')
set(gca,'XTickLabel','0|5|10|15|20|25|30')
subplot(3,1,2);
plot(Y(:,2),'b');
xlabel('t');
ylabel('x_{2}')
set(gca,'XTickLabel','0|5|10|15|20|25|30')
subplot(3,1,3);
plot(Y(:,3),'b');
xlabel('t');
ylabel('x_{3}')
set(gca,'XTickLabel','0|5|10|15|20|25|30')

figure
subplot(2,2,1);
plot(Y(:,1),Y(:,2),'b');
xlabel('x_{1}');
ylabel('x_{2}');
subplot(2,2,2);
plot(Y(:,1),Y(:,3),'b');
xlabel('x_{1}');
ylabel('x_{3}');
subplot(2,2,3);
plot(Y(:,2),Y(:,3),'b');
xlabel('x_{2}');
ylabel('x_{3}');
 
figure
plot(Y(:,1),Y(:,2),'b');
xlabel('x_{1}');
ylabel('x_{2}');
figure
plot(Y(:,1),Y(:,3),'b');
xlabel('x_{1}');
ylabel('x_{3}');
figure
plot(Y(:,2),Y(:,3),'b');
xlabel('x_{2}');
ylabel('x_{3}');

figure
plot3(Y(:,1),Y(:,2),Y(:,3),'b');
xlabel('x_{1}');
ylabel('x_{2}');
zlabel('x_{3}');
% T
 
%
