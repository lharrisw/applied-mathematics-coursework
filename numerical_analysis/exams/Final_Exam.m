%% Exercise 4
clc;
close all;
clear;

f = @(x,y) cos(10*(x.^2 + y.^2));
g = @(s) cos(10*s);

I = integral2(f,-1,1,-1,1);
Q = quad_gauss_exam(g,50,eps('double'));
E = I - Q;

%% Exercise 5
clc;
close all;
clear;

rhs = @(t,y) y.^3 - 5*y.^2 + 6*y + cos(t); % rhs of ODE
t0 = 0; % initial time
tf = 20; % final time

% step-sizes
n = linspace(0,6,7); 
h = 1./(2.^n);

E = zeros(1,7);

% initial condition
y0 = 0;

% plots
figure
hold on
subplot(4,2,[1 2])
[y,t] = ode_solver(rhs,t0,tf,h(1),y0); % Backward Euler solution
[tt,yy] = ode45(rhs,[t0 tf],y0); % ode45 solution
plot(t,y,'b-'); % plot of my solution
plot(tt,yy,'r.'); % plot of ode45 solution
xlabel('$t$','interpreter','latex');
ylabel('$y$','interpreter','latex');

E(1) = abs(norm(y - yy,'inf'));
% plots with stepsizes h = 0.5 to 0.015625
for k = 0:5
    [y,t] = ode_solver(rhs,t0,tf,h(k+2),y0);
    [tt,yy] = ode45(rhs,[t0 tf],y0);
    subplot(4,2,k+3)
    E(k+2) = abs(norm(y - yy,'inf'));
    hold on
    plot(t,y,'b-');
    plot(tt,yy,'r.');
    xlabel('$t$','interpreter','latex');
    ylabel('$y$','interpreter','latex');
end

n = linspace(0,15,16);
h = 1./(2.^n);
E = zeros(1,16);

for k = 1:16
    [y,t] = ode_solver(rhs,t0,tf,h(k),y0);
    [tt,yy] = ode45(rhs,[t0 tf],y0);
    E(k) = abs(norm(y - yy,'inf'));
end

x = abs(log(h(3:end)));
y = abs(log(E(3:end)));
p = polyfit(x,y,1);
xx = linspace(0,11,100);
pp = polyval(p,xx);

figure
hold on
grid on
plot(xx,pp,'b-');
plot(x,y,'bo');
legend('Fit','Log-Log Data');
title(sprintf("Order = %1.5f",p(1)));

% ODE
function [y,t] = ode_solver(rhs,t0,tf,h,y0)
    t = t0:h:tf; % time vector
    y = zeros(size(t)); % solution vector
    y(1) = y0; % initial condition
    
    for k = 1:length(t)-1
        ynew = y(k) + rhs(t(k),y(k))*h; % forward euler
        
        % Using forward euler's and newton's methods to approximate y(n+1)
        % with 100 iterations
        for j = 1:100
            g = ynew - y(k) + h*(ynew^3 -5*ynew^2 + 6*ynew + cos(t(k+1)));
            dg = 1 + h*(3*ynew^2 - 10*ynew + 6);
            ynew = ynew - g/dg;
        end
        
        % computing the the sequence
        y(k+1) = y(k) + h*rhs(t(k+1),ynew);
    end
end