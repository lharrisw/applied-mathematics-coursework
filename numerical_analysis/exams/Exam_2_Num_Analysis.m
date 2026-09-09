%% Exercise 1 Exam 2
clc
clear

N = 1e3;
m = 10000;
t = zeros(size(1,N));
tt = t;
err = zeros(m,N);
tol = 1e-4;

for k = 1:N
    x = abs(randn(m,1));
    tic
    y = 1./sqrt(x);
    t(k) = toc;
    tic
    yy = InvSqrt(x);
    tt(k) = toc;
    err(:,k) = abs(y - yy)';
end

fails = 0;

for j = 1:N
    if err(j) > tol
        fails = fails + 1;
    end
end

fail_rate = fails/N * 100;

mean_t = mean(t);
std_t = std(t);

mean_tt = mean(tt);
std_tt = std(tt);

mean_err = mean(err);
std_err = std(err);

vals = [y yy];

vals(1:5,:);
ratio = mean_tt/mean_t;

%% Exercise 1 Exam 2 Proof of quadratic convergence

y = @(s) cos(s) - s.^3;
dy = @(s) -sin(s) - 3*s.^2;
y1 = @(s) cos(2*s) - sin(2*s)./2 + s.^4;
dy1 = @(s) -2*sin(2*s) - cos(2*s) + 4*s.^3;

x = newt_method(y,dy,0.5,eps('double'));
yy = y(x);
err = abs(x - fzero(y,0.5))';
x1 = newt_method(y1,dy1,0.4,eps('double'));
err1 = abs(x1 - fzero(y1,0.4))';
yy1 = y1(x1)';

errs = [err(1:end-1) err1];


%% Exercise 2 Exam 2

y = @(x) LegendreP(100,x);
f = @(x) y(x).*cos(500*x);

I = integral(f,-1,1);
Q = quad_gauss(f,1000,1e-4);

% t = ones(size(1,10000));
% tt = t;
% 
% 
% 
% for k = 1:10000
%     tic
%     I = integral(f,-1,1);
%     t(k) = toc;
%     tic
%     Q = quad_gauss(f,10000,1e-4);
%     tt(k) = toc;
% end

% tic
% I = integral(f,-1,1);
% toc
% tic
% Q = quad_gauss(f,1000,1e-4);
% toc
% 
% mean_t = mean(t);
% std_t = std(t);
% 
% mean_tt = mean(tt);
% std_tt = std(tt);

%% Exercise 3 Exam 2
close;
clear;
clc;

% Solution

p0 = @(x) LegendreP(0,x);
p1 = @(x) LegendreP(1,x);
p2 = @(x) LegendreP(2,x);
p3 = @(x) LegendreP(3,x);

% Computed from Maple
c = [(exp(1)^2 - 1)/(2*exp(1)) 3/exp(1) (5*exp(1)^2 - 35)/(2*exp(1))...
    (-35*exp(1)^2 + 259)/(2*exp(1))];

% Computed from my function
[d,G] = coeffs(4);

x = linspace(-1,1,20);
f = @(x) c(1)*p0(x) + c(2)*p1(x) + c(3)*p2(x) + c(4)*p3(x);
g = @(x) d(1)*p0(x) + d(2)*p1(x) + d(3)*p2(x) + d(4)*p3(x);

figure
hold on
plot(x,exp(x),'b-');
plot(x,f(x),'ro');
plot(x,g(x),'k*');
legend('$y=e^{x}$','Maple Interpolant','My Interpolant','interpreter','latex');
xlabel('$x$','interpreter','latex');
ylabel('$y$','interpreter','latex');
title('Interpolating the Exponential');
set(gca,'fontsize',20);


%% Exercise 4 Exam 2 Part 2

N = 1000;
n = 49;
x = linspace(1,2,N);

g = @(x) 1;
h = zeros(size(x));

for k = 0:n
    h = h + integral(@(x) ex2_fun(x,k).*g(x),1,2).*ex2_fun(x,k);
end

figure
hold on
plot(x,h);
plot(x,ones(size(x)));
legend('$h_{49}$','y = g(x)','interpreter','latex');
xlabel('$x$','interpreter','latex');
ylabel('$y$','interpreter','latex');
title('Interpolation of $y = g(x)$','interpreter','latex');
set(gca,'fontsize',20);

%% Exercise 4 Exam 2 Part 3

N = 7;
x = linspace(0,1,10000);
u = 1./x;
q = zeros(N+1,10000);

for k = 0:N
    q(k+1,:) = (1/sqrt(2*k+1)).*u.*ex2_fun(u,k);
    hold on
    plot(x,q(k+1,:));
end