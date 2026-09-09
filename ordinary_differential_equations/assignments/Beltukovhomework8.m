%% Exercise 1 Part a

tt = linspace(0,12,300);
dydt = @(t,y) [y(2); -4*y(2) - 3*y(1)];
y0 = [0;1];
[t,y] = ode45(dydt,[0 12],y0);
x = (-0.5)*exp(-3*tt) + (0.5)*exp(-tt);

figure
hold on
plot(t,y(:,1),'ro');
plot(tt,x,'b-');
xlabel('time (sec)');
ylabel('x(t)');
title('Exercise 1 Part a')
legend('Numerical Solution','Symbolic Solution');

%% Exercise 1 Part b

t = linspace(0,12,300);
dydt = @(t,y) [y(2); -4*y(2) - 4*y(1)];
y0 = [0;1];
[t,y] = ode45(dydt,[0 12],y0);
x = t.*exp(-2*t);

figure
hold on
plot(t,y(:,1),'ro');
plot(t,x,'b-');
xlabel('time (sec)');
ylabel('x(t)');
title('Exercise 1 Part b')
legend('Numerical Solution','Symbolic Solution');

%% Exercise 1 Part c

t = linspace(0,3*pi);
dydt = @(t,y) [y(2); -4*y(2) - 10*y(1)];
y0 = [0;1];
[t,y] = ode45(dydt,[0 3*pi],y0);
x = (sqrt(6)/6).*exp(-2*t).*sin(sqrt(6)*t);

figure
hold on
plot(t,y(:,1),'ro');
plot(t,x,'b-');
xlabel('time (sec)');
ylabel('x(t)');
title('Exercise 1 Part c')
legend('Numerical Solution','Symbolic Solution');

%% Exercise 2

t = linspace(0,6*pi);
a = 1;
b = 1;
dydt = @(t,y) [y(2);a*cos(t) + b*sin(t)];
[t,y] = ode45(dydt,[0 6*pi],[0;0]);
x = a + b.*t - a*cos(t) - b*sin(t);
x1 = (a + b*t);

figure
hold on
plot(t,y(:,1),'ro');
plot(t,x,'b-');
xlabel('time (sec)');
legend('Numerical Solution','Symbolic Solution');
ylabel('x(t)');
ylim([0 .5]);
title('Exercise 2')

%% Exercise 3

t = linspace(0,10*pi,200);
x0 = 1;
x1 = .5;
x2 = 3;
x3 = .1;
dydt = @(t,y) [y(2);y(3);y(4);-y(3)];
y0 = [x0;x1;x2;x3];
[t,y] = ode45(dydt,[0 10*pi],y0);
x = -x2*cos(t) + -x3*sin(t) + (x1 + x3)*t + (x0 + x2);

figure
hold on
plot(t,y(:,1),'ro');
plot(t,x,'b-');
ylim([0 20]);
legend('Numerical Solution','Symbolic Solution');
title('Exercise 3')