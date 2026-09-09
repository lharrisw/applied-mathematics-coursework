%% Exercise 1

r = linspace(0,10,1000);

% plotting the kernel

a = 2;
b = 5;
c = 7;
x = 3;

% standard kernel plot

figure
plot(r,quad_kern(x,r,a,b,c));
xlabel('$r$','interpreter','latex')
ylabel(sprintf('$K(%1.f,r)$',x),'interpreter','latex');
title('Kernel of Remainder: Quadratic Interpolation','Interpreter',"latex");
set(gca,'fontsize',20);

% plot of first derivative

r1 = r(:,1:end-1);

figure
plot(r1,diff(quad_kern(x,r,a,b,c)));
xlabel('$r$','interpreter','latex')
ylabel(sprintf('$K''(%1.f,r)$',x),'interpreter','latex');
title('First Derivative of the Kernel','Interpreter','latex');
set(gca,'fontsize',20);

% plot of the second derivative

r2 = r1(:,1:end-1);

figure
plot(r2,diff(diff(quad_kern(x,r,a,b,c))));
xlabel('$r$','interpreter','latex')
ylabel(sprintf('$K''''(%1.f,r)$',x),'interpreter','latex');
title('Second Derivative of the Kernel','Interpreter',"latex");
set(gca,'fontsize',20);

%% Exercise 2

N = 5; % number of nodes and degree of interpolating polynomial
D = N-1; % degree of interpolant polynomial

x = linspace(0,1,N); % interval from 0 to 1 of equidistant nodes
y = exp(x); % function in question
p = polyfit(x,y,D); % coefficients for interpolant polynomial of degree D
xx = linspace(0,1,100); % finer grid for standard plot
yy = exp(xx); % function evaluated on finer grid

% Graphing

figure
hold on
plot(xx,polyval(p,xx),'r--'); % plot of the lagrange interpolating polynomial
plot(xx,yy,'k'); % plot of the function
plot(x,y,'o'); % plot of nodes
xlabel('$x$','Interpreter','latex'); % x-axis
ylabel('$y$','Interpreter','latex'); % y-axis
legend('Interpolant Polynomial','Function','Nodes');
title(sprintf('Lagrange Polynomial with Degree = %1.f, N = %1.f',D,N),'interpreter','latex');
set(gca,'fontsize',23);

for k = 1:9
    subplot(3,3,k);
    hold on
    
    u = linspace(0,1,k+1); % interval from 0 to 1 of equidistant nodes
    v = exp(u); % function in question
    pp = polyfit(x,y,k); % coefficients for interpolant polynomial of degree D
    uu = linspace(0,1,100); % finer grid for standard plot
    vv = exp(xx); % function evaluated on finer grid
    plot(uu,polyval(pp,uu),'r--'); % plot of the lagrange interpolating polynomial
    plot(uu,vv,'k'); % plot of the function
    plot(u,v,'o'); % plot of nodes
    xlabel('$x$','Interpreter','latex'); % x-axis
    ylabel('$y$','Interpreter','latex'); % y-axis
end

%% Exercise 3

N = 11; % number of nodes and degree of interpolating polynomial
D = N-1; % degree of interpolant polynomial

x = linspace(0,1,N); % interval form 0 to 1 of equidistant nodes
y = humps(x); % function in question
p = polyfit(x,y,D); % coefficients for interpolant polynomial of degree D
xx = linspace(0,1,100); % finer grid for standard plot
yy = humps(xx); % function evaluated on finer grid

% Graph

figure
hold on
plot(xx,polyval(p,xx),'r--'); % plot of the lagrange interpolating polynomial
plot(xx,yy,'k'); % plot of the function
plot(x,y,'o'); % plot of nodes
xlabel('$x$','Interpreter','latex'); % x-axis
ylabel('$y$','Interpreter','latex'); % y-axis
ylim([min(y),max(y)]);
legend('Interpolant Polynomial','Function','Nodes');
title(sprintf('Lagrange Polynomial with Degree = %1.f, N = %1.f',D,N),'interpreter','latex');
set(gca,'fontsize',23);

for k = 1:9
    subplot(3,3,k);
    hold on
    
    u = linspace(0,1,k+1); % interval from 0 to 1 of equidistant nodes
    v = humps(u); % function in question
    pp = polyfit(x,y,k); % coefficients for interpolant polynomial of degree D
    uu = linspace(0,1,100); % finer grid for standard plot
    vv = humps(xx); % function evaluated on finer grid
    plot(uu,polyval(pp,uu),'r--'); % plot of the lagrange interpolating polynomial
    plot(uu,vv,'k'); % plot of the function
    plot(u,v,'o'); % plot of nodes
    xlabel('$x$','Interpreter','latex'); % x-axis
    ylabel('$y$','Interpreter','latex'); % y-axis
end
