%% Exercise 1: Solving the Heat Equation w/ General Periodic Boundary Conditions

a = 2;
L = 2*a;
K = 1;

% Semidiscretization

N = 80; 
h = L/(N+1);
x = linspace(-a,a,N);

% Forming the matrix with Periodic Boundary Conditions

A = toeplitz([-2 1 zeros(1,N-2)]);
A(1,N) = 1;
A(N,1) = 1;

% Eigenvalue decomposition for A

[V,D] = eig(A); 
lambda = K*diag(D)/(h^2);

% Solution construction

f = @init_temp;
u0 = f(x)';
t = 0.1;
C = V\u0;
u = V*(exp(lambda*t).*C);
tt = linspace(0,100,N);

% Analytic Solution

en = zeros(N,1);
proj = zeros(length(x),N);

for n = 1:N
    g = @(x) cos(pi.*x.*n/a);
    gg = @(x) f(x).*g(x);
    proj(:,n) = (1/a)*quadgk(gg,-1,1)*g(x);
    en(n) = exp(-((n^2 * pi^2)/(a^2))*t);
end

uu = (1/(2*a))*quadgk(f,-1,1) + proj*en;

% Plotting the Analytic and Numerical Solutions

figure
hold on
plot(x,u0,'r-'); % Initial Condition
p = plot(x,u,'b-'); % Numerical Solution
plot(x,uu,'ko'); % Analytic Solution
title(sprintf('The Solution at t = %1.4f and N = %1.4f',[t N]));
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
%         pause(0.001);
%     end
% end

%% Exercise 2

l = linspace(-5000,0,100000);
y1 = @(l) l.*sin(sqrt(-l)); 
y2 = @(l) sqrt(-l).*cos(sqrt(-l));
dif = @(l) y2(l)-y1(l);  % Difference between both graphs

% Determining the eigenvalues

x = -1000000:20:0;
n = numel(x); % Produces the number of elements in x
z = zeros(1,n); % Creating an array

for j = 1:numel(x)
    z(j) = fzero(dif,x(j)); % Determines the eigenvalue from the dif
end

z = uniquetol(z,1e-6); % Provides the unique eigenvalues
z = sort(abs(z)); % Sorts the eigenvlues

% Plotting the first 9 eigenfunctions

x1 = linspace(0,10,1000);

for k = 1:9
    subplot(3,3,k);
    hold on
    plot(x1,(cos(x1*sqrt(z(k))) + sin(x1*sqrt(z(k))))); % The eigenfunctions
    title(sprintf('lambda = %1.4f',z(k)));
end

% Determining the rate of growth of eigenvectors

x2 = linspace(0,320);
p = polyfit([1:length(z)],z,2);
pp = polyval(p,x2);

% Confirming that the rate of growth is quadratic 

dz = diff(z); 
p1 = polyfit([1:length(z)-1],dz,1); 
pp1 = polyval(p1,x2(1:end-1));

yy = log(z);
pfit = polyfit(log([2:length(z)]),yy(2:end),1);
pval = polyval(pfit,log([2:length(z)]));

figure
hold on
plot(log([2:length(z)]),yy(2:end),'b-');
plot(log([2:length(z)]),pval,'ro');
title(sprintf('p = [%1.4f %1.4f]',[pfit]));

figure
hold on
plot(z,'b-');
plot(x2,pp,'ro');

figure
hold on
plot(dz,'b-');
plot(x2(1:end-1),pp1,'ro');


% figure
% hold on
% plot(l,y1(l),'r-');
% plot(l,y2(l),'b-');
% 
% figure
% hold on
% plot(l,diff(l));
% yline(0);

%% Exercise 3

% Crappy Code

f = @(c) ((2.*(21.*c(1).^2 + 70.*c(1).*c(3) + 35.*c(2).^2 + 105.*c(3).^2 + 15 - 42.*c(2)))/(105)).^2;

guess = [10000 20000 30000];
min = fminsearch(f,guess);

for g = 1:10000
    guess(g,:) = [g g*2 g*3];
    min(g,:) = fminsearch(f,guess(g,:));
end

min = mean(min);

x = linspace(-1,1);

p = @(x) x.^3;
q = @(c,x) c(1).*x.^2 + c(2).*x + c(3);

q = @(x) q(min,x);

h = @(x) p(x)-q(x);

norm = sqrt(quadgk(h,-1,1));

%% Exercise 3

% Better Code

% basis of V

v1 = @(x) 1;
v2 = @(x) x;
v3 = @(x) x.^2;

f = @(x) x.^3;

% Gram Matrix

g = zeros(3,3);

g(1,1) = integral(@(x) v1(x).*v1(x),-1,1,'arrayvalued',true);
g(1,2) = integral(@(x) v1(x).*v2(x),-1,1);
g(1,3) = integral(@(x) v1(x).*v3(x),-1,1);
g(2,1) = quadgk(@(x) v2(x).*v1(x),-1,1);
g(2,2) = quadgk(@(x) v2(x).*v2(x),-1,1);
g(2,3) = quadgk(@(x) v2(x).*v3(x),-1,1);
g(3,1) = quadgk(@(x) v3(x).*v1(x),-1,1);
g(3,2) = quadgk(@(x) v3(x).*v2(x),-1,1);
g(3,3) = quadgk(@(x) v3(x).*v3(x),-1,1);

% Solution Vector

b = zeros(3,1);

b(1,1) = quadgk(@(x) v1(x).*f(x),-1,1);
b(2,1) = quadgk(@(x) v2(x).*f(x),-1,1);
b(3,1) = quadgk(@(x) v3(x).*f(x),-1,1);

% Solution
c = g\b;