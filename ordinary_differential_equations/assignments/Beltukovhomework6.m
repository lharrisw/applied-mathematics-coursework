%% ODE Investigation
clear
close all
clc

t = linspace(0,2*pi,200);
dxdt = @(t,x) cos(3*t) - 2*x;
x0 = 10;
T = 2*pi;
[t,x] = ode45(dxdt,[0 T],x0);
ind1 = t > 3;
t1 = t(ind1);
x1 = x(ind1);
A1 = [cos(3*t1) sin(3*t1)];
b1 = linsolve(A1,x1);
x2 = A1*b1;
ind2 = t < 3;
t2 = t(ind2);
x3 = x(ind2);
A2 = [cos(3*t2) sin(3*t2)];
b2 = x3 - A2*b1;
c = log(b2);
p = polyfit(t2,c,1);
x4 = b1(1)*cos(3*t) + b1(2)*sin(3*t) + exp(p(2))*exp(p(1)*t);
x5 = (128/13)*exp(-2*t) + (2/13)*cos(3*t) + (3/13)*sin(3*t);

figure
hold on
plot(t,x,'ro');
plot(t,x4,'b-');
plot(t,x5,'*');
xlabel('t');
ylabel('x');
legend('Numerical Solution','Investigative Solution','Symbolic Solution');

figure
hold on
plot(t1,x1,'b-');
plot(t1,x2,'ro');
xlabel('t > 3');
ylabel('x(t>3)');
title('The steady part of the solution');

figure
hold on
plot(t2,c,'ro');
plot(t2,polyval(p,t2),'b-');
xlabel('t');
ylabel('log(t<3)');
title('The transient part of the solution');

%% ODE Investigation 2
close all
clc
clear

t = linspace(0,10,200);
dxdt = @(t,x) t - 3*x;
x0 = -5:1:5;
T = 10;
[t,x] = ode45(dxdt,[0 T],x0);
%ind1 = t > 3;
%t1 = t(ind1);
%x1 = x(ind1);
%A1 = [cos(3*t1) sin(3*t1)];
%b1 = linsolve(A1,x1);
%x2 = A1*b;
%ind2 = t < 3;
%t2 = t(ind2);
%x3 = x(ind2);
%A2 = [cos(3*t2) sin(3*t2)];
%b2 = x3 - A2*b1;
%c = log(b2);
%p = polyfit(t2,c,1);
%x4 = b1(1)*cos(3*t) + b1(2)*sin(3*t) + exp(p(2))*exp(p(1)*t);
%x5 = (128/13)*exp(-2*t) + (2/13)*cos(3*t) + (3/13)*sin(3*t);

figure
hold on
plot(t,x,'b-');
%plot(t,x4,'b-');
%plot(t,x5,'*');
xlabel('t');
ylabel('x');
legend('Solutions to equation (3)');
%legend('Numerical Solution','Investigative Solution','Symbolic Solution');

figure
hold on
plot(t1,x1,'b-');
plot(t1,x2,'ro');
xlabel('t > 3');
ylabel('x(t>3)');
title('The steady part of the solution');

figure
hold on
plot(t2,c,'ro');
plot(t2,polyval(p,t2),'b-');
xlabel('t');
ylabel('log(t<3)');
title('The transient part of the solution');