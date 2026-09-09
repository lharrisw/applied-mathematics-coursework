t = linspace(0,3*pi);
dxdt = @(t,x) [x(2) + 1 - cos(x(1)); -4*x(1) - sin(x(2))];
x0 = [0;1];
[t,x] = ode45(dxdt,[0 3*pi],x0);
p = [1 1 4];
l = roots(p);
A = [1 1;l(1) l(2)];
b = [0;1];
C = linsolve(A,b);
x1 = C(1)*exp(l(1)*t) + C(2)*exp(l(2)*t);

figure
hold on
plot(t,x(:,1),'b-');
plot(t,x1,'ro');