%% Exercise 1
clear
clc
close all

N = 5000;

one_norm = zeros(1,N);

for n = 1:N
    dir_ker = @(s) abs(sin((n+0.5)*s)./sin(s/2));
    one_norm(n) = integral(dir_ker,0,2*pi)/(2*pi);
end

x = log(1:N);
y = one_norm;

p = polyfit(x,y,1);
pp = polyval(p,x);

% for k = 1:4
%     subplot(2,2,k)
%     hold on
%     plot(s,(sin((k^3+0.5)*s)./sin(s/2)));
%     title(sprintf('$$ N = %1.0f $$',[k^3]),'interpreter','latex');
% end

figure
hold on,
plot(x,y,'b-');
plot(x,pp,'ro');
xlabel('$$log(N)$$','interpreter','latex');
ylabel('$$L^{1}$$-Norm of $$D_{N}$$','interpreter','latex');
legend('Data','Fit');
title(sprintf('Constant = %f',p(1)));
set(gca,'FontSize',18);

%% Exercise 2

%f = @(x) abs(x);
f = @(x) abs(x).^2; % C1 function
%f = @(x) abs(x).^3; % C2 function
%f = @(x) abs(x).^4; % C3 function
%f = @(x) abs(x).^500; % C499 function
%f = @(x) exp(x);
%f = @(x) exp(-1i*x);
N = 50;
nn = 50;

% c = zeros(1,N);
% for n = 1:N
%     c(n) = integral(@(x)exp(-1i*n*x).*f(x),-pi,pi)/(2*pi);
% end

for k = 1:620
    g = @(x) sin(x);
    for l = 1:nn
        c1(l) = integral(@(x)exp(-1i*l*x).*g(x),-pi,pi)/(2*pi);
    end
    uu = 1:620;
    vv = abs(c1);
    ind = vv>eps('single');
    uu = log(uu(ind));
    vv = log(vv(ind));
    p1(k,:) = polyfit(uu,vv,1);
end

% logkk = log(1:620)';
% logp1 = log(p1(:,1));
% pk = polyfit(logkk(end-400:end),logp1(end-400:end),1);
% ppk = polyval(pk,logkk(end-400:end));
% 
% 
% figure
% hold on
% plot(logkk,logp1,'b-');
% plot(logkk(end-400:end),ppk(end-400:end),'ro');
% ylabel('Log of Rate of Decay');
% xlabel('Log of Differentiability');
% legend('Decay vs. Differentiability','Data Fit');
% title(sprintf('Order of Convergence: p = %1.4f, N = %1.0f',[pk(1),nn]));
% set(gca,'FontSize',18);

figure
hold on
plot((1:620),abs(p1(:,1)),'b-');
xlabel('k');
ylabel('Rate of Decay');
title(sprintf('The Rate of Decay r(k,N=%1.0f)',nn));
set(gca,'FontSize',18);

% xx = 1:N;
% yy = abs(c);
% ind = yy>eps('single');
% xx = log(xx(ind));
% yy = log(yy(ind));
% p = polyfit(xx(end-20:end),yy(end-20:end),1);
% 
% figure
% hold on
% plot(xx(end-20:end),yy(end-20:end),'b-');
% plot(xx(end-20:end),polyval(p,xx(end-20:end)),'ro');
% xlabel('log(N)');
% ylabel('Log of Fourier Coefficients');
% title(sprintf('Rate of Decay of |x|^{2} = [%1.4f %1.4f]',p));
% set(gca,'FontSize',18);


%% Exercises 3 and 4

%f = @(x) abs(x);
%f = @(x) abs(x).^2;
%f = @(x) abs(x).^3;
%f = @(x) abs(x).^4;
%f = @(x) sign(x);
%f = @(x) sin(x);
%f = @(x) 1.*(x<0.5) + 0.*(x>0.5);
N = 50;

x = linspace(0,1,100);
sum_b = zeros(size(x));
%xx = zeros(N,N);

for l = 1:N
    for k = 1:l
        for j = 0:k
            b = @(x) f(j/N).*nchoosek(N,j).*(x.^j).*((1-x).^(N-j));
        end
    end
    sum_b = sum_b + b(x);
    R = (sum_b - f(x))/l;
    E(k) = max(abs(R));
end

xx = log(1:N-1);
yy = log(E(1:end-1));
p = polyfit(xx,yy,1);

% Plot of the Polynomial againss the Actual Function

figure
hold on
plot(x,f(x),'b-');
plot(x,sum_b,'ro');
xlabel('x');
ylabel('f(x)');
legend('Graph','Bernstein Fit');
title('Bernstein Fit of f(x) = sin(x)');
set(gca,'FontSize',18);

% The Order of Convergence

figure
hold on
plot(xx,yy,'b-');
plot(xx,polyval(p,xx),'ro');
xlabel('log(N)');
ylabel('log(E)');
legend('Error','Fit');
title(sprintf('Order of Covergence = [%1.3f %1.3f]',p));
set(gca,'FontSize',18);

% Plot of the Error
figure
hold on
plot(exp(xx),exp(yy),'b-');
xlabel('N');
ylabel('Error');
title('Plot of the Error');
set(gca,'FontSize',18);

%% Exercise 5

L = pi;
N = 100;
x = linspace(0,2*L,N);

f = @(x) sign(pi - x);

err_l2 = zeros(size(1:N));
err_inf = zeros(size(1:N));
ss = zeros(length(x),length(x));
SS = zeros(length(x),length(x));

% Infinity-norm and L2-norm convergence tests

for n = 1:N
    ss = zeros(size(x));
    for m = 1:n
        g = @(x) cos(m*x);
        g1 = @(x) sin(m*x);
        gg = @(x) g(x).*f(x);
        gg1 = @(x) g1(x).*f(x);
        ss = ss + (1/L)*quadgk(gg,0,2*L)*g(x) + (1/L)*quadgk(gg1,0,2*L)*g1(x);
    end
    ss = ss + quadgk(f,0,2*L)*ones(size(ss))./(4*L);
    R = abs(f(x) - ss);
    err_l2(n) = norm((R.^2).*(2*L)/n,2);
    err_inf(n) = max(abs(R));
end

% L2-norm fit

logN = log(1:N);
logE_l2 = log(err_l2);
p_l2 = polyfit(logN(end-20:end),logE_l2(end-20:end),1);
order_l2 = p_l2(1);

% Infinity-norm "fit"

logN = log(1:N);
logE_inf = log(err_inf);
p_inf = polyfit(logN(end-20:end),logE_inf(end-20:end),1);
order_inf = p_inf(1);

% Plot of L2-norm convergence

figure
hold on
plot(logN,logE_l2,'ro');
plot(logN,polyval(p_l2,logN),'b-');
xlabel('Log of N');
ylabel('Log of Error');
title(sprintf('$$L^{2}$$-Norm Convergence: %f',order_l2),'interpreter','latex');
set(gca,'FontSize',18);

% Plot of Infinity-norm convergence

figure
hold on
plot(logN,logE_inf,'ro');
plot(logN,polyval(p_inf,logN),'b-');
xlabel('Log of N');
ylabel('Log of Error');
title(sprintf('Infinity Norm Convergence: %f',order_inf));
set(gca,'FontSize',18);

% Plot of function and Fourier Approximation

figure
hold on
plot(x,f(x),'b-');
plot(x,ss,'ro');
xlabel('x');
ylabel('f(x)');
legend('Function','Fourier Approximation');
title('The Function and Fourier Approximation');
set(gca,'FontSize',18);

%% Exercise 6 Legendre Polynomials

N = 500;
x = linspace(-1,1,1000);

f = @(x) sign(x);

s = zeros(size(x));
ss = zeros(size(x));

for n = 1:N
    g = @(x) f(x).*legendre_n(x,(n-1));
    proj = ((2*(n-1)+1)/2)*quadgk(g,-1,1).*legendre_n(x,(n-1));
    s = s + proj;
    R = f(x) - s;
    E(n) = max(abs(R));
end

% for j = 1:6
%     subplot(2,3,j);
%     ss = zeros(size(x));
%     for nn = 1:(2^j)
%         gg = @(x) f(x).*legendre_n(x,(nn-1));
%         proj1 = ((2*(nn-1)+1)/2)*quadgk(gg,-1,1).*legendre_n(x,(nn-1));
%         ss = ss + proj1;
%     end
%     hold on
%     plot(x,f(x),'b-');
%     plot(x,ss,'r-');
%     title(sprintf('N = %1.0f',nn));
%     set(gca,'FontSize',18);
% end

% xx = 1:N;
% yy = E;
% p = polyfit(xx,yy,1);
% pp = polyval(p,xx);
% 
% % Fit
% 
% figure
% hold on
% plot(xx,yy,'b-');
% plot(xx,pp,'ro');
% xlabel('N');
% ylabel('Error');
% title(sprintf('Legendre Polynomial Convergence: %1.3f',p(1)));
% set(gca,'FontSize',18);

figure
hold on
plot(x,f(x));
plot(x,s);