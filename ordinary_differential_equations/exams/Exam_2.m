%% Exercise 1
t = linspace(0,4*pi);
dxdt = @(t,x) [x(2); cos(t) + cos(3*t) + cos(5*t) - 4*x(1)];
x0 = [0;0];
T = 4*pi;
[t,x] = ode45(dxdt,[0 T],x0);
xp = 0.1667*(exp(i*1*t) + exp(-i*1*t)) - 0.1000*(exp(i*3*t) + exp(-i*3*t)) - 0.0238*(exp(i*5*t) + exp(-i*5*t));
A = ones(2,2);
A(2,:) = [2*i -2*i];
b = [-0.0429; 0];
C = linsolve(A,b);
xc = C(1)*exp(2*i*t) + C(2)*exp(-2*i*t);
x1 = xc + xp;

figure
hold on
plot(t,x(:,1),'ro');
plot(t,x1,'b-');
xlabel('t');
ylabel('x');
legend('Numerical Solution','Symbolic Solution');
title('Exercise 1 Part c');

%% Exercise 4

w = 2;
t = linspace(0,5*pi);
dxdt = @(t,x) [x(2);sin(w*t) - 4*x(1)];
x0 = [0;0];
T = 4*pi;
[t,x] = ode45(dxdt,[0 T],x0); % numerical solution
p = (i*w)^2 + 4;
x1 = sin(w*t)/p - w*sin(2*t)/(2*p); % for omega = 1
x2 = -0.5*(0.5*t.*cos(2*t) - 0.25*sin(2*t)); % for omega = 2

figure
hold on
plot(t,x(:,1),'ro'); % for omega = 1
plot(t,x1,'b-');
xlabel('t');
ylabel('x');
legend('Numerical Solution','Symbolic Solution');

figure
hold on
plot(t,x(:,1),'ro'); % for omega = 2
plot(t,x2,'b');
xlabel('t');
ylabel('x');
legend('Numerical Solution','Symbolic Solution');