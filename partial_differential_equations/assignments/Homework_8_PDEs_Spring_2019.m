%% Inverse Problem Heat Equation

% Forming a Symmetric Grid

L = pi;
N = 100; 
h = L/N;
xx = -L:h:L;

% Heat Equation Solution on Grid

t = 1;
%f = @(x) 3*x;
%f = @(x) 100*cos(3*x) - 90*sin(2*x);
%f = @(x) exp(-x.^2);
%f = @sinc;
%f = @humps;
%f = @sign;
%f = @sawtooth;
%f = @(x) cos(4*x + t);

u = zeros(size(xx)); 

for k = 1:length(xx)
    hh = @(y) exp(-((xx(k)-y).^2)/(4*t))/sqrt(4*pi*t); % Representer
    u(k) = quadgk(@(y) hh(y).*f(y),-L,L); % Coefficients (solution at specified time)
end

% Gram matrix

gram = zeros(length(xx),length(xx));

for ii = 1:length(xx)
    for jj = 1:length(xx)
        gram(ii,jj) = quadgk(@(y) (1/(4*pi*t))*exp(-((xx(ii)-y).^2)/(4*t))...
            .*exp((-(xx(jj)-y).^2)/(4*t)),-L,L);
    end
end

% Minimum-norm solution

reps = zeros(length(xx),length(xx));

for n = 1:length(xx)
    reps(:,n) = exp((-(xx(n)-xx).^2)/(4*t))/sqrt(4*pi*t); % Matrix of representers as columns
end

% Solving an Ill conditioned system

c = linsolve(gram,u');
c = gram\u';
ff = reps*c;

% Utilizing the Psuedoinverse

c1 = pinv(gram)*u';
ff1 = reps*c1;

figure
hold on
plot(xx,ff,'b-');
plot(xx,f(xx),'ro');

figure
hold on
plot(xx,ff1,'b-');
plot(xx,f(xx),'ro');

% xlabel('$x$','interpreter','latex');
% ylabel('$u(x,0)$','interpreter','latex');
% legend('Minimum-Norm Solution','Solution with Pinv','Actual Initial Condition','interpreter','latex');
% title(sprintf('Initial Condition Approximation, N = %1.0f',N),'interpreter','latex');
% set(gca,'fontsize',18);

%% Inverse problem Wave Equation

% Forming a Symmetric Grid

L = pi;
N = 300; 
h = (2*L)/N;
xx = -L:h:L;

% Heat Equation Solution on Grid

f = @(x) 3*x;
%f = @(x) 100*cos(3*x) - 90*sin(2*x);
%f = @(x) exp(-x.^2); % Initial Data
%f = @sinc;
%f = @humps;
%f = @sign;
f = @sawtooth;
t = 1;

u = zeros(size(xx)); 

for k = 1:length(t)
    u(k) = quad2d(@(z,y) exp(1i*z*(x(k)-z)).*cos(z*t).*f(y),-L,L,-L,L)/(2*pi)...
        + quad2d(@(z,y) exp(1i*z*(x(k)-z)).*sin(z*t).*f(y)./z,-L,L,-L,L)/(2*pi);
end

% Gram matrix

gram = zeros(length(xx),length(xx));

for ii = 1:length(xx)
    for jj = 1:length(xx)
        gram(ii,jj) = quadgk(@(z) exp(1i*z*(xx(ii)-y)).*cos(z*t)/(2*pi)...
            .*exp(1i*z*(x(jj)-y)).*cos(z*t),-L,L)/(2*pi);
    end
end

% Minimum-norm solution

reps = zeros(length(xx),length(xx));

for n = 1:length(xx)
    reps(:,n) = exp((-(xx(n)-xx).^2)/(4*t))/sqrt(4*pi*t); % Matrix of representers as columns
end

c = linsolve(gram,u'); 

ff = reps*c;

figure
hold on
plot(xx,ff,'b-');
plot(xx,f(xx),'ro');
xlabel('$x$','interpreter','latex');
ylabel('$u(x,0)$','interpreter','latex');
legend('Minimum-Norm Solution','Actual Initial Condition','interpreter','latex');
title(sprintf('Initial Condition Approximation: $f = sawtooth(x)$, N = %1.0f',N),'interpreter','latex');
set(gca,'fontsize',18);