%% Exercise 1

% Function evaluations

u = [0,1,2];
v = [-10,12,8];

% Lagrange interpolation

xx = linspace(0,2);
p = polyfit(u,v,2);
pp = polyval(p,xx);

% Functions

f = @(x) p(1)*x.^2 + p(2)*x + p(3); % interpolated quadratic
dfdx = @(x) 2*p(1)*x + p(2);

f1 = @(x) sin(x); % my function
df1dx = @(x) cos(x);

% Roots

[x1,x2] = newt_method(f,dfdx,0.25,5);
root1 = fzero(f,0.25);

[xx1,xx2] = newt_method(f1,df1dx,3*pi/4,5);
root2 = fzero(f1,3*pi/4);

% Graph of interpolant

figure
hold on
grid on 
plot(xx,pp);
plot(u,v,'bo');
plot(x1,f(x1),'b*',x2,f(x2),'r*');
xlabel('$x$','interpreter','latex');
ylabel('$y$','interpreter','latex');
legend('Interpolant','Interpolation Points','$x_{1}$','$x_{2}$','interpreter','latex');
title(sprintf('Root of $f(x)$: $x_{1} = %1.4f$, $x_{2} = %1.4f$',[x1,x2]),'interpreter','latex');
set(gca,'fontsize',20);

% Graph of my function

figure
hold on
grid on
plot(0:0.01:2*pi,f1(0:0.01:2*pi));
plot(xx1,f1(xx1),'b*',xx2,f1(xx2),'r*');
xlabel('$x$','interpreter','latex');
ylabel('$y$','interpreter','latex');
legend('$f(x) = \sin(x)$','$x_{1}$','$x_{2}$','interpreter','latex');
title(sprintf('Root of $f(x) =$ sin$(x)$: $x_{1} = %1.4f$, $x_{2} = %1.4f$',[xx1,xx2]),'interpreter','latex');
set(gca,'fontsize',20);

%% Exercise 4

t = linspace(-1,1);
K = @(t) (((1-t).^4)/4 - ((-1/sqrt(3) - t).^3).*(-1 <= t <= -1/sqrt(3))...
    - ((1/sqrt(3) - t).^3).*(-1 <= t <= 1/sqrt(3)))/6; 

plot(t,K(t));

%% Exercise 5 Part 3

n = 500;
a = -1;
b = 1;
E = zeros(size(1:n));

c = randn(1,11);
%f = @(x) polyval(c,x);
f = @(x) x.^4;
I = integral(f,a,b);
hh = zeros(size(1:n));
x = linspace(-1,1,n);

for k = 1:n
    h = (b-a)/(k); % multiply k by 2 for comp_simp and by 3 for comp_simp38
    %Q = comp_simp38(f,a,b,3*k);
    %Q = trap_quad(f,a,b,k);
    %Q = mid_quad(f,a,b,k);
    Q = comp_simp(f,a,b,2*k);
    E(k) = abs(I-Q);
    hh(k) = h;
end

u = log(hh);
v = log(E);
p = polyfit(u,v,1);
txt1 = 'Midpoint Rule';
txt2 = 'Trapezoid Rule';
txt3 = 'Composite Simpson''s Rule';
txt4 = 'Composite Simpson''s 3/8 Rule';

figure
hold on 
plot(u,v,'ro');
plot(u,polyval(p,u));
xlabel('$\log(h)$','interpreter','latex');
ylabel('$\log(E)$','interpreter','latex');
legend('Data','Fit');
title(sprintf('The order of the %s: m = %1.4f',txt4,p(1)),'interpreter','latex');
set(gca,'fontsize',20);

figure
hold on
plot(x,f(x))
title('A random polynomial of degree 10.','interpreter','latex');
set(gca,'fontsize',20);

%% Exercise 5 part 5

n = 5000;
a = -1;
b = 1;

%fail_asimp38 = 0;
fail_asimp13 = 0;
fail_quad = 0;

for j = 1:n
    %tol = 1/(1000*n);
    tol = 1e-14;
    c = randn(1,11);
    I = polyval(polyint(c),b) - polyval(polyint(c),a);
    f = @(x) polyval(c,x);
    %Q1 = a_simp38(f,a,b,tol);
    Q1 = a_simp13(f,a,b,tol);
    Q2 = quad(f,a,b,tol);
    E1 = abs(I-Q1);
    E2 = abs(I-Q2);
    if E1 > tol
        %fail_asimp38 = fail_asimp38 + 1;
        fail_asimp13 = fail_asimp13 + 1;
    end
    if E2 > tol
        fail_quad = fail_quad + 1;
    end
end

%rate1 = fail_asimp38/n;
rate1 = fail_asimp13/n;
rate2 = fail_quad/n;


%% Exercise 5 Part 6

x = [1e-1 1e-2 1e-3 1e-4 1e-5 1e-6 1e-7 1e-8 1e-9 1e-10];
y = [43 61 115 157 223 331 463 769 955 1591];

logx = log(x);
logy = log(y);
[p,S] = polyfit(logx,logy,1);
pp = polyval(p,logx);

figure
hold on
plot(logx,logy,'bo');
plot(logx,pp,'b-');
xlabel('$\log(tol)$','interpreter','latex');
ylabel('$\log(evals)$','interpreter','latex');
title(['Log of Function' ...
    ' Evaluations vs. Log of Tolerance'],'interpreter','latex');
legend('Data','Fit','interpreter','latex');
set(gca,'fontsize',20);