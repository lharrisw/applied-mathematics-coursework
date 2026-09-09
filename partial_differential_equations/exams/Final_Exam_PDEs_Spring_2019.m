%% Exercise 1

N = 50;
m = 5;
n = 25;

a = 2;
%a = 4;
%a = 8;
%a = 16;

b = 1;

[r,theta] = meshgrid(linspace(0,a,N),linspace(0,2*pi,N));

% Right-hand side

f = @(r,t) 1.*(r<b) + 0.*(r>b);

% Computing first n zeros of derivatives of desired bessel eigenfunctions

z = zeros(n+1,n);

for j = 0:n
    for k = 1:n
        z(j+1,k) = BessDerivZerosBisect2(j,k,eps('double')); 
    end
end

z(1,:) = circshift(z(1,:),1);
z(1,1) = 0;

% Solution for Neumann condition

u = zeros(size(r));

% for ii = 1:n-1
%     ll = (z(1,ii)/a).^2;
%     aa0 = (integral2(@(rr,tt) f(rr,tt).*besselj(0,sqrt(ll).*rr).*rr,0,a,0,2*pi)...
%             ./integral(@(rr) ((besselj(0,sqrt(ll).*rr)).^2).*rr,0,a))/(2*pi);
%     u = u + (-1)*besselj(0,sqrt(ll)*r).*aa0;
% end
% 
% for mm = 1:m
%     for nn = 1:n
%         l = (z(mm+1,nn)/a).^2;
%         aa = (integral2(@(rr,tt) f(rr,tt).*besselj(mm,sqrt(l).*rr).*cos(mm*tt).*rr,0,a,0,2*pi)...
%             ./integral(@(rr) ((besselj(mm,sqrt(l).*rr)).^2).*rr,0,a))/(pi);
%         bb = (integral2(@(rr,tt) f(rr,tt).*besselj(mm,sqrt(l).*rr).*sin(mm*tt).*rr,0,a,0,2*pi)...
%             ./integral(@(rr) ((besselj(mm,sqrt(l).*rr)).^2).*rr,0,a))/(pi);
%         u = u + (-1)*besselj(mm,sqrt(l)*r).*(aa*cos(mm*theta) + bb*sin(mm*theta));
%     end
% end

% Verification of Initial Condition

ff = zeros(size(r));

for ii = 1:n-1
    ll = (z(1,ii)/a).^2;
    aa0 = (integral2(@(rr,tt) f(rr,tt).*besselj(0,sqrt(ll).*rr).*rr,0,a,0,2*pi)...
            ./integral(@(rr) ((besselj(0,sqrt(ll).*rr)).^2).*rr,0,a))/(2*pi);
    ff = ff + besselj(0,sqrt(ll)*r).*aa0;
end

for mm = 1:m
    for nn = 1:n
        l = (z(mm+1,nn)/a).^2;
        aa = (integral2(@(rr,tt) f(rr,tt).*besselj(mm,sqrt(l).*rr).*cos(mm*tt).*rr,0,a,0,2*pi)...
            ./integral(@(rr) ((besselj(mm,sqrt(l).*rr)).^2).*rr,0,a))/(pi);
        bb = (integral2(@(rr,tt) f(rr,tt).*besselj(mm,sqrt(l).*rr).*sin(mm*tt).*rr,0,a,0,2*pi)...
            ./integral(@(rr) ((besselj(mm,sqrt(l).*rr)).^2).*rr,0,a))/(pi);
        ff = ff + besselj(mm,sqrt(l)*r).*(aa*cos(mm*theta) + bb*sin(mm*theta));
    end
end

% Plot Solutions

% figure
% mesh(r.*cos(theta),r.*sin(theta),u); 
% xlabel('$x$','Interpreter','latex');
% ylabel('$y$','Interpreter','latex');
% zlabel('$u(r,\theta)$','Interpreter','latex');
% title(sprintf('Poisson Equation Solution: $a$ = %1.0f',a),'interpreter','latex');
% set(gca,'fontsize',18);

figure
subplot(1,2,1)
mesh(r.*cos(theta),r.*sin(theta),f(r,theta)); 
xlabel('$x$','Interpreter','latex');
ylabel('$y$','Interpreter','latex');
zlabel('$f$','Interpreter','latex');
title('Actual Initial Condition');
set(gca,'fontsize',18);
zlim([-0.2 1.2]);

subplot(1,2,2)
mesh(r.*cos(theta),r.*sin(theta),ff); 
xlabel('$x$','Interpreter','latex');
ylabel('$y$','Interpreter','latex');
zlabel('$f$','Interpreter','latex');
title('Eigenfunction Initial Condition');
set(gca,'fontsize',18);
zlim([-0.2 1.2]);

%% PDEs Final Exam Question 1 Redone Correctly

clear;
clc;
close all;

% Outer disk radius
% a = 2;
a = 4;
% a = 8;
% a = 16;

% inner disk radius
b = 1;

% formation of grid
N = 50;
[r,theta] = meshgrid(linspace(0,a,N),linspace(0,2*pi,N));

% Right-hand side 
f = @(r,t) 1.*(r<b)+0.*(r>=b);

% Computing first n zeros of derivatives of desired bessel eigenfunctions
n = 50; 
z = zeros(1,n);

for j = 1:n
    z(j) = BessDerivZerosBisect2(0,j,eps('single')); 
end

z(1,:) = circshift(z,1);
z(1,1) = 0;

% Projection of f onto L^2(Ω) and solution to new systems

constant = 0;
ff = zeros(size(r));
u = zeros(size(r))+constant;

for k = 2:n
    lambda = (z(k)/a)^2;
    c = integral2(@(rr,tt)f(rr,tt).*besselj(0,sqrt(lambda).*rr).*rr,0,b,0,2*pi)...
        /integral2(@(rr,tt)besselj(0,sqrt(lambda).*rr).^2.*rr,0,a,0,2*pi);
    A = -c/lambda;
    ff = ff + c*besselj(0,sqrt(lambda)*r); % projection of f onto L^2(Ω)
    u = u + A*besselj(0,sqrt(lambda)*r); % solution
end

% Plot of new solutions

figure
mesh(r.*cos(theta),r.*sin(theta),u,'AlphaData',u,...
    'FaceAlpha','flat')
% mesh(r.*cos(theta),r.*sin(theta),u); 
xlabel('$x$','Interpreter','latex');
ylabel('$y$','Interpreter','latex');
zlabel('$u(r,\theta)$','Interpreter','latex');
title('Solution to $\nabla u = \tilde{f} =$ proj$_{\Omega}(f)$'...
      ,sprintf('$a$ = %1.0f',a),'interpreter','latex');
set(gca,'fontsize',18);
axis tight

% Plot of f and proj(f)

figure
subplot(1,2,1)
mesh(r.*cos(theta),r.*sin(theta),f(r,theta)); 
xlabel('$x$','Interpreter','latex');
ylabel('$y$','Interpreter','latex');
zlabel('$f$','Interpreter','latex');
title('The Right Hand Side, $f$','interpreter','latex');
set(gca,'fontsize',18);
axis tight
 
subplot(1,2,2)
mesh(r.*cos(theta),r.*sin(theta),ff); 
xlabel('$x$','Interpreter','latex');
ylabel('$y$','Interpreter','latex');
zlabel('proj$_{\Omega}(f)$','Interpreter','latex');
title('Projection of $f$ onto $\Omega$','interpreter','latex');
set(gca,'fontsize',18);
axis tight

%% Exercise 2

t = [1 2 3 4]; % Desired times

x = linspace(-5,5,100); % Length of finite domain

u = zeros(length(x),1); % Solution vectors
uu = u;

% Initial condition

f = @(x) exp(-2*x.^2);

% Transform Solution

figure
for tt = 1:length(t)
    for j = 1:length(x)
        u(j,:) = integral2(@(w,y) exp(1i*w.*(x(j)-y+t(tt))).*f(y),-10,10,-10,10)/(2*pi);
    end
    hold on
    subplot(2,2,tt)
    plot(x,u,'b-');
    xlabel('Rod Position','interpreter','latex');
    ylabel('$u(x,t)$','interpreter','latex');
    legend(sprintf('$t$ = %1.0f',t(tt)),'interpreter','latex');
    title('Transform Solution','interpreter','latex');
    set(gca,'fontsize',18);
    ylim([-0.05 1.05]);
end

% Explicit Solution

figure
for j = 1:length(t)
    uu = exp(-2*(x+t(j)).^2); 
    hold on
    subplot(2,2,j);
    plot(x,uu,'r-');
    xlabel('Rod Position','interpreter','latex');
    ylabel('$u(x,t)$','interpreter','latex');
    legend(sprintf('$t$ = %1.0f',t(j)),'interpreter','latex');
    title('Explicit Solution','interpreter','latex');
    set(gca,'fontsize',18);
    ylim([-0.05 1.05]);
end

%% Exercise 3

N = [50 100 200 400 800 1600]; 

for j = 1:length(N)
    
    % Semi-discretization
    
    L = pi;
    h = L/(N(j)+1); % step-size
    x = (1:N(j))*h; % 1D grid
    
    % Generating the matrix
    
    f = (-4*cos(4*x)./(1+x)); % function evaluated on grid
    a = 1./(1+x); 
    d = eye(N(j)+1).*((16*h^2) - 2); % main diagonal
    da = diag(1-a*h,1); % subdiagonal above
    db = diag(1+a*h,-1); % subdiagonal below
    A = d + da + db; % Full matrix with Dirichlet conditions
    A(:,N(j)+1) = []; % Delete extra column
    A(N(j)+1,:) = []; % Delete extra row
    
    % Solving the linear system
    
    u = A\(f*h^2)'; 
    
    % Plotting to show convergence
    
    hold on
    plot(x,u);
    xlabel('$x$','interpreter','latex');
    ylabel('$u(x)$','interpreter','latex');
    title('Solution');
    set(gca,'fontsize',18);
end

legend('N = 50','N = 100','N = 200','N = 400','N = 800','N = 1600');