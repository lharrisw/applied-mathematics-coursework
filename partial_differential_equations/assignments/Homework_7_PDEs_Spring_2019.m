%% Poisson Kernel and Solution to Laplace Equation on a Disk

% Generating the grid

N = 100;
[r,theta] = meshgrid(linspace(0,1,N),linspace(0,2*pi,N));

% Initial Conditions

%g = @(p) log(p).*p;
g = @(p) cos(4*p) + sin(2*p);

% Solution

u = zeros(size(r));

for n = 0:N
    u = u + (1/pi)*integral(@(t) g(t).*cos(n*t),0,2*pi).*cos(n*theta).*r.^(abs(n))...
        + (1/pi)*integral(@(t) g(t).*sin(n*t),0,2*pi).*sin(n*theta).*r.^(abs(n)); 
end

% Plotting in Cartesian Coordinates

figure
mesh(r.*cos(theta),r.*sin(theta),real(u)); 
xlabel('$x$','Interpreter',"latex");
ylabel('$y$','Interpreter',"latex");
zlabel('$u(x,y)$','Interpreter',"latex");
title('Laplacian on the Unit Disk');
set(gca,'fontsize',18);