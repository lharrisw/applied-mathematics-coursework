%% Exercise 1 part a
t = linspace(0,6*pi,200);
dxdt = @(t,x) cos(t) - 2*x;
x0 = 10;
T = 6*pi;
[tt,x] = ode45(dxdt,[0 T],x0);
x1 = (48/5)*exp(-2*t) + (2/5)*cos(t) + (1/5)*sin(t);

figure
hold on
plot(tt,x,'ro');
plot(t,x1,'b-');
xlabel('time');
ylabel('x');
legend('Numerical Solution','Symbolic Solution');
title('Exercise 1 Part (a)');
%% Exercise 1 part b
t = linspace(0,5,200);
dxdt = @(t,x) exp(-t) - 2*x;
x0 = 10;
T = 5;
[t,x] = ode45(dxdt,[0 T],x0);
x1 = 9*exp(-2*t) + exp(-t);

figure
hold on
plot(t,x,'ro');
plot(t,x1,'b-');
xlabel('time');
ylabel('x');
legend('Numerical Solution','Symbolic Solution');
title('Exercise 1 Part (b)');
%% Exercise 1 part c
t = linspace(0,6*pi,200);
dxdt = @(t,x) 2*cos(t) + 3*exp(-t) - 2*x;
x0 = 10;
T = 6*pi;
[t,x] = ode45(dxdt,[0 T],x0);
x1 = (31/5)*exp(-2*t) + 3*exp(-t) + (4/5)*cos(t) + (2/5)*sin(t);

figure
hold on
plot(t,x,'ro');
plot(t,x1,'b-');
xlabel('time');
ylabel('x');
legend('Numerical Solution','Symbolic Solution');
title('Exercise 1 Part (c)');
%% Exercise 2
w = 2;
m = 1.5;
r = 3;
A = 4;
t = linspace(0,4*pi,200);
T = 4*pi;
x10 = 3;
x20 = 8;
f = @(t,x) [x(2); (A*sin(w*t) - r*x(2))/m];
[t,x] = ode45(f,[0 T],[x10;x20]);
C = (x20 + (A*m*w)/(((m^2)*(w^2)) + (r^2)));
D = x10 + (m/r)*C + (A*r)/(((m^2)*(w^3)) + ((r^2)*(w)));
x1 = D + (-m/r)*C*exp((-r/m)*t) - (A*r*cos(w*t))/(((m^2)*(w^3)) + ((r^2)*(w))) - (A*m*sin(w*t))/(((m^2)*(w^2)) + (r^2));


figure
hold on
plot(t,x(:,1),'ro');
plot(t,x1,'b');
xlabel('time');
ylabel('x');
legend('Numerical Solution','Symbolic Solution');
title('Exercise 2');
ylim([0 10]);

%% Exercise 3

A = [1 2;2 4];
b = [3;6];