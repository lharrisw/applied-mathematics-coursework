%% Exercise  3
clear
clc
close all

L = 1;
n = 6;
x = linspace(0,L,n);
h = x(2)-x(1);
f = @(t) cos(pi*t/h);

for i = 1:n
    for j = 1:n
        area(i,j) = quadgk(@(t) f(t),x(i),x(j));
    end
end

%% Exercise 4

h = 0.01; % Generating known data
x = -1:h:1; 
x0 = 10;
[m,n] = size(x);
hh = linspace(-1,1,n);


y = x.*exp(x); % The function for the approximations
yy = @(t) t.*exp(t); % The function for the error analysis


d2ydx2 = @(x) 2*exp(x) + x.*exp(x); % The true derivative of y


% Central Differences (CD)
for l = 2:n-1
    dy(l) = (y(l+1) - 2*y(l) + y(l-1))/(h^2);
end


dyy = (yy(x0 + hh) - 2*yy(x0) + yy(x0 - hh))./(hh.^2); % Error approx for CD
errc = abs(dyy - d2ydx2(x0));
logh = log(hh(hh>0)); % Use positive values of hh due to log
logerrc = log(errc);
p = polyfit(logh,logerrc(hh>0),1); % Use same values of hh to ensure same vector length
pfit = polyval(p,logh);


% Forward Differences (FD)
for j = 1:n-2
    dy1(j) = (y(j+2) - 2*y(j+1) + y(j))/(h^2);
end


dyy1 = (yy(x0 + 2*hh) - 2*yy(x0 + hh) + yy(x0))./(hh.^2); % Error approx for FD
errf = abs(dyy1 - d2ydx2(x0));
logerrf = log(errf);
pp = polyfit(logh,logerrf(hh>0),1); 
ppfit = polyval(pp,logh);


figure
hold on 
plot(x,d2ydx2(x),'ro'); % Plots the actual derivative 
plot(x(1:end-2),dy1,'b-'); % Plots the forward differences approx.
plot(x(1:end-1),dy,'k-'); % Plots the central differences approx.
xlabel('x');
ylabel('y');
legend('True Derivative','Forward Differences Approximation',...
    'Central Differences Approximation');
title('Approximations of $\displaystyle\frac{d^2y}{dx^2} = 2e^x + xe^x;~h = 0.01$',...
    'interpreter','latex');
hold off

figure
hold on
plot(logh,logerrc(hh>0),'ro'); % Plot of stepsize vs. CD error; Maintain same vector length
plot(logh,pfit,'b-'); % Plot of the fit
xlabel('Log of Stepsize');
ylabel('Log of Error');
title(sprintf('Central Differences; p = [%f %f], h = 0.01',p));
legend('Error Data','Fit');

figure
hold on
plot(logh,logerrf(hh>0),'ro'); % Plot of stepsize vs. FD error; Maintain same vector length
plot(logh,ppfit,'b-'); % Plot of the fit
xlabel('Log of Stepsize');
ylabel('Log of Error');
title(sprintf('Foward Differences; p = [%f %f], h = 0.01',pp));
legend('Error Data','Fit');

%% Exercise 5: Solving the heat equation, discrete case
clear
clc
close all

L = 1;
K = 1;
 
% Semidiscretization

N = 1000; 
h = L/(N+1);
x = (1:N)*h;

% Forming the matrix with Dirichlet Conditions

A = toeplitz([-2 1 zeros(1,N-2)]);
A(1,1) = -1;
A(N,N) = -1;

% Eigenval decomp for A

[V,D] = eig(A); 
lambda = K*diag(D)/(h^2);

% Solution construction

f = @humps;
u0 = f(x)';
t = 3;
C = V\u0;
u = V*(exp(lambda*t).*C);
tt = linspace(0,1,N);

% Plotting the solution

figure
hold on
p = plot(x,u,'b-');
plot(x,u0,'r-');
xlabel('x');
ylabel(sprintf('u(x,%1.4f)',t));
title(sprintf('Solution at time t = %f',t));

% Configuring the eigenvectors for plotting

% [lambda,l] = sort(abs(diag(D)));
% V = V(:,l);
% lambda = lambda./(h^2);

% Animating the solution

for t = tt
    u = V*(exp(lambda*t).*C);
    if ishandle(p)
        set(p,'ydata',u);
        title(sprintf('t=%1.4f',t),'FontSize',20);
        pause(0.1);
    end
end

% for k = 1:9
%     subplot(3,3,k);
%     plot(x,V(:,k)./V(1:k),'r.');
%     hold on
%     ylim([-1 1]);
% end