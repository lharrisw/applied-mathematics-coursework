%% Exercise 2

a = randn;
b = randn;

t = linspace(a,b,100);

K = @(t) ((((a-t).^3).*(a+2*b-3*t))/12).*(a <= t & t < (a+b)/2)...
    + ((((b-t).^3).*(2*a+b-3*t))/12).*((a+b)/2 <= t &  t <= b);

figure
plot(t,K(t));
xlabel('$t$','interpreter','latex');
ylabel('$K(t)$','interpreter','latex');
title('Plot of the Kernel','interpreter','latex');
set(gca,'fontsize',20);



%% Exercise 8 and 10

clear;
clc;

n = 1000;
a = -pi;
b = pi;
E = zeros(size(1:n));

%f = @(x) 1./(2+abs(sin(x)));
%f = @(x) 1./(2+sin(x));
%f = @(x) abs(asec(x))./(1+abs(cos(x))); % 
%f = @(x) x.^4;
%f = @(x) sin(x);
%f = @humps;
%f = @(x) (sin(x)).^(1/10);
%f = @(x) abs((sin(x)).^(1/10));
f = @(x) x.^(1/5);
I = integral(f,a,b); 
hh = zeros(size(1:n));
x = linspace(a,b,n);
pp = zeros(size(1:n));

for k = 1:n
    h = (b-a)/(k); % multiply k by 2 for comp_simp and by 3 for comp_simp38
    %Q = comp_simp(f,a,b,2*k);
    Q = comp_trap(f,a,b,k);
    E(k) = abs(I-Q);
    hh(k) = h;
end

% y = @(b,x) b(1).*exp(-b(2).*x) + b(3).*exp(-b(4).*x) + b(5);
% xx  = u;
% yx = v;
% OLS = @(b) sum((y(b,xx) - yx).^2);           % Ordinary Least Squares cost function
% opts = optimset('MaxFunEvals',50000, 'MaxIter',10000);
% B = fminsearch(OLS,[300; 0.028; 100; 0.0047; 4], opts);       % Use 1fminsearch0 to minimise the 0OLS0 function
% fcnfit = y(B,xx);                            % Calculate function with estimated parameters
% figure
% hold on
% plot(xx, yx, '*b');
% plot(xx, fcnfit, '-r');

u = log(hh);
% u = log(1./(1:n));
v = log(E);
p = polyfit(u,v,1);
%txt = 'Simpson''s 1/3 Rule';
txt = 'Trapezoid Rule';

figure
hold on 
plot(u,v,'ro');
plot(u,polyval(p,u));
xlabel('$\log(h)$','interpreter','latex');
ylabel('$\log(E)$','interpreter','latex');
legend('Data','Fit');
title(sprintf('The order of the %s: m = %1.4f',txt,p(1)),'interpreter','latex');
set(gca,'fontsize',20);

figure
hold on
plot(x,f(x));
%title('Plot of $1/(2+|\sin(x)|)$','interpreter','latex');
%title('Plot of $1/(2+\sin(x))$','interpreter','latex');
%title('Plot of humps$(x)$','interpreter','latex');
xlabel('$x$','interpreter','latex');
ylabel('$y$','interpreter','latex');
set(gca,'fontsize',20);

%% Exercise 9

n = 10000;
a = 0;
b = 1;

fail_trap = 0;
fail_mid = 0;
k = abs(randn);

tic
for j = 1:n
    tol = 10^(-5*k);
    c = randn(1,6);
    I = polyval(polyint(c),b) - polyval(polyint(c),a);
    f = @(x) polyval(c,x);
    Q1 = comp_trap(f,a,b,j);
    Q2 = comp_mid(f,a,b,j);
    E1 = abs(I-Q1);
    E2 = abs(I-Q2);
    if E1 > tol
        fail_trap = fail_trap + 1;
    end
    if E2 > tol
        fail_mid = fail_mid + 1;
    end
end
toc

rate1 = fail_trap/n;
rate2 = fail_mid/n;
ratio = rate1/rate2;


%% Extra stuffs

[x,y] = meshgrid(linspace(-1,1,100));

f = @(x,y) (pi^2)/2 - 4*atanh(x)./x - 4*atanh(y)./y + 1./(1-x.^2) ...
    + 1./(1-y.^2) + 2./(1-(x.*y));

f1 = @(x) (pi^2)/2 - 4*atanh(x)./x + 4*atanh(-x)./x + 1./(1-x.^2) ...
    + 1./(1-x.^2) + 2./(1+x.^2);

minval1 = fminsearch(f1,0.6);


s = linspace(-0.75,0.75,1000);
g = @(s) -atanh(s);
h = @(s) (s.^7 + 3*s.^5 - s.^3 -s)./((s.^4 - 1).^2);
fun = @(s) h(s)-g(s);
miny = fminsearch(fun,-0.5);

minval = fminsearch(@(v) f(v(1),v(2)),[0.6 0.6]);

figure
mesh(x,y,f(x,y));
xlabel('$x$','interpreter','latex');
ylabel('$y$','interpreter','latex');
zlabel('$f(x,y)$','interpreter','latex');

figure
hold on
plot(s,g(s));
plot(s,h(s));