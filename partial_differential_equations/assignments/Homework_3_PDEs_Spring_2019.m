%% Exercise 1: Solving the Heat Equation w/ Periodic Boundary Conditions\
clear
clc
close all

L = 2*pi;
K = 1;

% Semidiscretization

N = 300; 
h = L/(N+1);
x = linspace(-pi,pi,N);

% Forming the matrix with Periodic Boundary Conditions

A = toeplitz([-2 1 zeros(1,N-2)]);
A(1,N) = 1;
A(N,1) = 1;

% Eigenvalue decomposition for A

[V,D] = eig(A); 
lambda = K*diag(D)/(h^2);

% Solution construction

%f = @humps;
%f = @(x) sin(3*x);
%f = @(x) abs(x);
f = @sign; 
u0 = f(x)';
t = 1.9;
C = V\u0;
u = V*(exp(lambda*t).*C);
tt = linspace(0,3,N);

% Analytic Solution

for n = 1:N
    g = @(x) sin(x.*n);
    gg = @(x) f(x).*g(x);
    g1 = @(x) cos(x.*n);
    gg1 = @(x) f(x).*g1(x);
    proj(:,n) = (1/pi)*quadgk(gg,-pi,pi)*g(x);
    proj1(:,n) = (1/pi)*quadgk(gg1,-pi,pi)*g1(x);
    en(n) = exp(-(n^2)*t);
end

uu = (1/(2*pi))*quadgk(f,-pi,pi) + (proj + proj1)*en';

% Plotting the Analytic and Numerical Solutions

figure
hold on
plot(x,u0,'r-'); % Initial Condition
p = plot(x,u,'b-'); % Numerical Solution
%plot(x,uu,'ko'); % Analytic Solution
title(sprintf('The Solution at t = %1.4f',t));
legend('Initial Condition','Numerical Solution','Analytic Solution');
ylabel('u(x,t)');
xlabel('x');
set(gca,'FontSize',18);

% Animating the solution

% for t = tt
%     u = V*(exp(lambda*t).*C);
%     if ishandle(p)
%         set(p,'ydata',u);
%         title(sprintf('t=%1.4f',t),'FontSize',20);
%         ylabel(sprintf('u(x,%1.4f)',t))
%         pause(0.000001);
%     end
% end


% for k = 1:9
%     subplot(3,3,k);
%     hold on
%     plot(x,sign(V(1,k)).*V(:,k),'r.');
%     plot(x,sqrt(2/(N+1))*sin(k*x));
% end

%% Exercise 2: Solving the Wave Equation w/ Periodic Boundary Conditions
clear
close all
clc

L = 2*pi;
K = 1;

% Semidiscretization

N = 200; 
h = L/(N+1);
x = linspace(-pi,pi,N);

% Forming the matrix with Periodic Boundary Conditions

A = toeplitz([-2 1 zeros(1,N-2)]);
A(1,N) = 1;
A(N,1) = 1;

% Eigenvalue decomposition for A

[V,D] = eig(A); 
lambda = K*diag(D)/(h^2);

% Solution construction

%f = @humps;
%f = @(x) sin(3*x);
f = @(x) abs(x);
%f = @sign;
u0 = f(x)';
t = 1;
C = V\u0;
u = V*(cos(sqrt(-lambda)*t).*C);
tt = linspace(0,10,N);

% Analytic Solution

for n = 1:N
    g = @(x) sin(x.*n);
    gg = @(x) f(x).*g(x);
    g1 = @(x) cos(x.*n);
    gg1 = @(x) f(x).*g1(x);
    proj(:,n) = (1/pi)*quadgk(gg,-pi,pi)*g(x);
    proj1(:,n) = (1/pi)*quadgk(gg1,-pi,pi)*g1(x);
    harm(n) = cos(n.*t);
end

uu = (1/(2*pi))*quadgk(f,-pi,pi) + (proj + proj1)*harm';

% Plotting the Analytic and Numerical Solutions

figure
hold on
plot(x,u0,'r-'); % Initial Condition
p = plot(x,u,'b-'); % Numerical Solution
plot(x,uu,'ko'); % Analytic Solution
title(sprintf('The Solution at t = %1.4f',t));
legend('Initial Condition','Numerical Solution','Analytic Solution');
ylabel('u(x,t)');
xlabel('x');
set(gca,'FontSize',18);

% Animating the solution

for t = tt
    u = V*(cos(sqrt(-lambda)*t).*C);
    if ishandle(p)
        set(p,'ydata',u);
        title(sprintf('t=%1.4f',t),'FontSize',20);
        ylabel(sprintf('u(x,%1.4f)',t))
        pause(0.000000000001);
    end
end