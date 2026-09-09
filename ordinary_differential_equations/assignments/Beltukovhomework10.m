%% Non-Homogeneous Matrix-Vector ODE

t = linspace(0,2*pi);
dxdt = @(t,x) [x(2);-x(1) - x(2) - sin(t)];
T = 2*pi;
x0 = [0;0];
[t,x] = ode45(dxdt,[0 T],x0);
A = [0 1;-1 -1];
I = [1 0;0 1];
[V,D] = eig(A);
lambda1 = D(1);
lambda2 = D(4);
a1b1 = (A^-1 + I)^-1 *[0;-1];
a2b2 = A*a1b1;
xp = a1b1*cos(t)' + a2b2*sin(t)';
B = ones(2,2);
B(2,:) = [lambda1 lambda2];
C = V^-1 * [1;0];
xc = C(1)*V(:,1)*exp(lambda1*t)' + C(2)*V(:,2)*exp(lambda2*t)';
x1 = xc + xp;

figure
hold on
plot(t,x(:,1),'ro');
plot(t,x(:,1),'b-');
xlabel('t');
ylabel('x');
legend('Numerical Solution','Matrix-Vector Solution');
title('Non-homogeneous Matrix-Vector ODE');