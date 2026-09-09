%% Exercise 1: Solving the Dirichlet Problem of the Heat Equation
clear
clc
close all

L = pi;
K = 1;

% Semidiscretization

N = 200; 
h = L/(N+1);
x = (1:N)*h;

% Forming the matrix with Dirichlet Conditions

A = toeplitz([-2 1 zeros(1,N-2)]);

% Eigenvalue decomposition for A

[V,D] = eig(A); 
lambda = K*diag(D)/(h^2);

% Solution construction

f = @humps;
%f = @sign;
u0 = f(x)';
t = 1;
C = V\u0;
u = V*(exp(lambda*t).*C);
tt = linspace(0,1,N);

% Configuring the eigenvectors for plotting

% [lambda,l] = sort(abs(diag(D)));
% V = V(:,l);
% lambda = lambda./(h^2);

% Analytic Solution

proj = zeros(N,N);
en = zeros(1,N);

for n = 1:N
    g = @(x) sin(pi.*x.*n/L);
    gg = @(x) f(x).*g(x);
    proj(:,n) = (2/L)*quadgk(gg,0,L)*g(x);
    en(n) = exp(-(n^2)*t);
end

uu = proj*en';

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

% Animating the solution

for t = tt
    u = V*(exp(lambda*t).*C);
    if ishandle(p)
        set(p,'ydata',u);
        title(sprintf('t=%1.4f',t),'FontSize',20);
        ylabel(sprintf('u(x,%1.4f)',t))
        pause(0.001);
     end
end


% for k = 1:9
%     subplot(3,3,k);
%     hold on
%     plot(x,sign(V(1,k)).*V(:,k),'r.');
%     plot(x,sqrt(2/(N+1))*sin(k*x));
%     legend('Eigenvectors from Eig','Eigenvectors from Derivation');
% end

%% Exercise 2 The Neumann Problem for the Wave Equation

clear;
close all;
clc;

L = pi;
K = 1;

% Semidiscretization

N = 500; 
h = L/(N+1);
x = (1:N)*h;

% Forming the matrix with Neumann Conditions

A = toeplitz([-2 1 zeros(1,N-2)]);

% Neumann conditions 
% A(1,1) = -1;
% A(N,N) = -1;

% Eigenvalue decomposition for A

[V,D] = eig(A); 
lambda = K*diag(D)/(h^2);

% Solution construction

f = @humps;
u0 = f(x)';
t = 1;
C = V\u0;
u = @(t) V*(cos(sqrt(-lambda)*t).*C);
tt = linspace(0,5,N);

% Analytic Solution

for n = 1:N
    g = @(x) cos(pi.*x.*n/L);
    gg = @(x) f(x).*g(x);
    proj(:,n) = (2/L)*quadgk(gg,0,L)*g(x);
    en(n) = cos(pi.*t.*n/L);
end

uu = (1/L)*quadgk(f,0,L) + proj*en';

% Plotting the solution

figure
hold on
% plot(x,u0,'r-'); % Initial Condition
p = plot(x,u(t),'b-'); % Analytic Solution
% plot(x,uu,'ko'); % Analytic Solution
title(sprintf('The Solution at t = %1.4f',t));
legend('Initial Condition','Numerical Solution','Analytic Solution');
ylabel('u(x,t)');
xlabel('x');

% for ii = 1:4
%     subplot(2,2,ii)
%     hold on
%     plot(x,u(exp(ii)^2 * 0.001),'b-');
%     plot(x,u0,'r-');
%     xlabel('x');
%     ylabel(sprintf('u(x,%1.4f)',exp(ii)^2 * 0.001));
%     title(sprintf('Solution at time t = %f',exp(ii)^2 * 0.001));
% end


% Configuring the eigenvectors for plotting

% [lambda,l] = sort(abs(diag(D)));
% V = V(:,l);
% lambda = lambda./(h^2);

% Animating the solution

for t = tt
    u = V*(cos(sqrt(-lambda)*t).*C);
    if ishandle(p)
        set(p,'ydata',u);
        title(sprintf('t=%1.4f',t),'FontSize',20);
        pause(0.001);
    end
end

% for k = 1:9
%     subplot(3,3,k);
%     hold on
%     plot(x,sign(V(1,k)).*V(:,k),'r.');
%     plot(x,sqrt(2/(N+1))*sin(k*x));
%     legend('Eigenvectors from Eig','Eigenvectors from Derivation');
% end

%% Fourier Error Analysis

% Conditions

L = 7*pi/3;
N = 1:100;
x = linspace(0,L,100);
w = pi/L;

% The functions

f = @sign;
%f = @(x) x;
%f = @(x) x.*(L - x);
%f = @(x) 1./(2 + cos(w*x));

err = zeros(size(N));
ss = zeros(length(x),length(x));

for n=N
    s = zeros(size(x));
    for ii=1:n
        g = @(x) cos(pi*ii/L*x);
        gf = @(x) g(x).*f(x);
        s = s + 2/L*quadgk(gf,0,L)*g(x);
    end
    s = s + 1/L*quadgk(f,0,L)*ones(size(s));
    SS(n,:) = s;
    R = f(x) - s;
    E(n) = max(abs(R));
end



logN = log(N);
logE = log(E);
p = polyfit(logN(end-20:end),logE(end-20:end),1);
order = p(1);

figure
hold on
plot(logN,logE,'ro')
plot(0:5,polyval(p,0:5),'b-')
xlabel('Log of N')
ylabel('Log of Error')
title(sprintf('Order of convergence: %f for f(x) = 1/(2 + cos(wx))',order));


%% Justin's Code for Exercise 3

L = 10;
%f = @(x) x;
%f = @(x) x.*(L-x);
%f = @(x) 1./(2+cos(pi/L*x));
f = @(x) sin(4*x).*exp(-0.5*x);
%f = @(x) sin(10*x);
x = linspace(0,L,100);

NN = 1:100;
E = zeros(size(NN));
SS = zeros(length(x), length(x));

for N=NN
    S = zeros(size(x));
    for n=1:N
        c = @(x) cos(pi*n/L*x);
        cf = @(x) c(x).*f(x);
        S = S + 2/L*quadgk(cf,0,L)*c(x);
    end
    S = S + 1/L*quadgk(f,0,L)*ones(size(S));
    SS(N,:) = S;
    R = f(x) - S;
    E(N) = max(abs(R));
end

logNN = log(NN);
logE = log(E);
p = polyfit(logNN(end-20:end),logE(end-20:end),1);
order = p(1);

 

figure
hold on
grid on
plot(logNN,logE,'ro')
plot(0:5,polyval(p,0:5),'b-')
xlabel('log(N)')
ylabel('log(R)')
title(sprintf('Order of convergence: %f', order))

% 

% figure
% N = 1:49:102;
% for i=1:3
%     subplot(3,1,i)
%     hold on
%     grid on
%     plot(x,f(x),'b-')
%     plot(x,SS(N(i),:),'ro');
%     ylim([min(f(x)) max(f(x))])
%     xlabel('x')
%     ylabel('f(x)')
%     title(sprintf('N = %d', N(i)))
% end

 

% figure
% hold on
% grid on
% plot(x,f(x),'b-')
% plt = plot(x,SS(1,:),'ro');
% ylim([min(f(x)) max(f(x))])
% xlabel('x')
% ylabel('f(x)')
% 
% for N=NN
%     s = SS(N,:);
%     if ishandle(plt)
%         set(plt,'ydata',s)
%         title(sprintf('N = %d', N))
%         pause(0.1)
%     end
% end

 

% figure
% subplot(2,1,1)
% hold on
% plot(x,f(x),'b-')
% plt = plot(x,SS(1,:),'ro');
% 
% for N=NN
%     subplot(2,1,1)
%     hold on
%     grid on
%     ylim([min(f(x)) max(f(x))])
%     s = SS(N,:);
%     if ishandle(plt)
%         set(plt,'ydata',s)
%         title(sprintf('N = %d', N))
%         pause(0.1)
%     end
%     subplot(2,1,2)
%     hold on
%     grid on
%     plot(logNN(N),logE(N),'ro')
%     xlim([0 5])
%     ylim([min(logE) max(logE)])
%     xlabel('log(N)')
%     ylabel('log(R)')
%     title(sprintf('Order of convergence: %f', order))
% end


% subplot(2,1,2)
% grid on
% hold on
% plot(0:5,polyval(p,0:5),'b-')