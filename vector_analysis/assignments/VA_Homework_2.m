%% Parameterize an Ellipse R2

t = linspace(0,2*pi,1000);
z = linspace(0,1,1000);
z1 = linspace(-1,0,1000);
w = linspace(0,sqrt(2),1000);
w1 = linspace(-sqrt(2),0,1000);
a = 1;
b = 1;
x0 = 0;
y0 = 0;
x = @(t) x0 + a*cos(t);
y = @(t) y0 + b*sin(t);

y1 = z-1;
y2 = -z1-1;
y3 = -z+1;
y4 = z1+1;

y5 = w - sqrt(2);
y6 = -w1 -sqrt(2);
y7 = -w + sqrt(2);
y8 = w1 + sqrt(2);

figure
hold on
grid on
plot(x(t),y(t));
plot(z,y1,'r');
plot(z1,y2,'r');
plot(z,y3,'r');
plot(z1,y4,'r');
plot(w,y5,'r');
plot(w1,y6,'r');
plot(w,y7,'r');
plot(w1,y8,'r');
axis equal

%% Parameterize a Sphere and Flux

% Parameterization
theta = linspace(0,2*pi);
phi = linspace(0,pi);
[t,p] = meshgrid(theta,phi); % create meshgrid of angles
r = 1;
x0 = 0; % x0 -> y0 make up the initial point
y0 = 0;
z0 = 0;
x = @(t,p) x0 + r*sin(p)*cos(t); % parameterization of sphere
y = @(t,p) y0 + r*sin(p)*sin(t);
z = @(t,p) z0 + r*cos(p);
norm = @(t,p) sqrt(x(t,p).^2 + y(t,p).^2 + z(t,p).^2);

% Electric field evaluated on the sphere
Ex = @(t,p) (1/(4*pi))*x(t,p)./((norm(t,p)).^3);
Ey = @(t,p) (1/(4*pi))*y(t,p)./((norm(t,p)).^3);
Ez = @(t,p) (1/(4*pi))*z(t,p)./((norm(t,p)).^3);
dydz = @(t,p) -(r^2)*cos(t)*sin(p)^2;
dxdz = @(t,p) (r^2)*sin(t)*sin(p)^2;
dxdy = @(t,p) -(r^2)*sin(p)*cos(p);
flux_form = @(t,p) Ex(t,p).*dydz(t,p) - Ey(t,p).*dxdz(t,p) + Ez(t,p).*dxdy(t,p);
flux = dblquad(flux_form,0,2*pi,0,pi);

surf(x(t,p),y(t,p),z(t,p));
colormap('spring');
shading interp
axis square

%% Parameterize Ellipsoid and Flux in R3

% Parameterization
theta = linspace(0,2*pi);
phi = linspace(0,pi);
[t,p] = meshgrid(theta,phi); % create meshgrid of angles
a = 2;
b = 6;
c = 3;
x0 = 1; % x0 -> y0 make up the initial point
y0 = 2;
z0 = 10;
x = @(t,p) x0 + a*sin(p)*cos(t); % parameterization of ellipsoid
y = @(t,p) y0 + b*sin(p)*sin(t);
z = @(t,p) z0 + c*cos(p);
norm = @(t,p) sqrt(x(t,p).^2 + y(t,p).^2 + z(t,p).^2);

% Electric field evaluated on the ellipsoid
Ex = @(t,p) (1/(4*pi))*x(t,p)./((norm(t,p)).^3);
Ey = @(t,p) (1/(4*pi))*y(t,p)./((norm(t,p)).^3);
Ez = @(t,p) (1/(4*pi))*z(t,p)./((norm(t,p)).^3);
dydz = @(t,p) -(b*c)*cos(t)*sin(p)^2;
dxdz = @(t,p) (a*c)*sin(t)*sin(p)^2;
dxdy = @(t,p) -(a*b)*sin(p)*cos(p);
flux_form = @(t,p) Ex(t,p).*dydz(t,p) - Ey(t,p).*dxdz(t,p) + Ez(t,p).*dxdy(t,p);
flux = dblquad(flux_form,0,2*pi,0,pi);

surf(x(t,p),y(t,p),z(t,p));
colormap('summer');
shading interp

%% Parameterizing a Torus

theta = linspace(0,2*pi,50);
psi = linspace(0,2*pi,50);
[p,t] = meshgrid(psi,theta); % create meshgrid of angles
r = 5; % radius from the center of the torus to the center of the tube
rho = 2; % radius of the center of the tube
x0 = 1; % x0 -> z0 comprise initial point
y0 = 1;
z0 = 3;
x = (r + rho.*cos(t)).*cos(p) + x0; % parameterization of torus
y = (r + rho.*cos(t)).*sin(p) + y0;
z = rho.*sin(t) + z0;


surf(x,y,z);
daspect([1 1 1]) % preserves the shape of torus
colormap('summer') % change color appearance 
shading interp
material metal
title('Torus')
xlabel('X');ylabel('Y');zlabel('Z');

%% Flux through a Torus

theta = linspace(0,2*pi);
psi = linspace(0,2*pi);
[p,t] = meshgrid(psi,theta); % create meshgrid of angles
r = 2; % radius from the center of the torus to the center of the tube
rho = 1; % radius of the center of the tube
x0 = 2; % x0 -> z0 comprise initial point
y0 = 0;
z0 = 0;
x = @(t,p) (r + rho.*cos(t)).*cos(p) + x0; % parameterization of torus
y = @(t,p) (r + rho.*cos(t)).*sin(p) + y0;
z = @(t,p) rho.*sin(t) + z0;
norm = @(t,p) sqrt(x(t,p).^2 + y(t,p).^2 + z(t,p).^2);

% Electric field evaluated on the ellipsoid
Ex = @(t,p) (1/(4*pi))*x(t,p)./((norm(t,p)).^3);
Ey = @(t,p) (1/(4*pi))*y(t,p)./((norm(t,p)).^3);
Ez = @(t,p) (1/(4*pi))*z(t,p)./((norm(t,p)).^3);
dydz = @(t,p) -(rho.*cos(t)).*(r.*cos(p) + rho.*cos(t).*cos(p));
dxdz = @(t,p) -(rho*cos(t)).*(-r.*sin(p) - rho.*cos(t).*cos(p));
dxdy = @(t,p) -rho.*sin(t).*cos(p).*(r.*cos(p) + rho.*cos(t).*cos(p)) ...
        + rho.*sin(t).*sin(p).*(-r.*sin(p) - rho.*cos(t).*sin(p));
flux_form = @(t,p) Ex(t,p).*dydz(t,p) - Ey(t,p).*dxdz(t,p) + Ez(t,p).*dxdy(t,p);
flux = dblquad(flux_form,0,2*pi,0,2*pi);

%% Exam 2 exercise 1

% The Field
h = 0.01;
n = 100;
m = 100;
[x,y] = meshgrid(-5:h:5);
for ii = 1:n
    u(ii,:) = ((x(ii,:)-ii).^2)/sqrt((x(ii,:)-ii).^2 + (y(ii,:)-ii).^2);
    for jj = 1:m
        v(jj,:) = ((y(jj,:)-jj).^2)/sqrt((x(jj,:)-jj).^2 + (y(jj,:)-jj).^2);
    end
end

% Divergence through the interior
divF = divergence(x,y,u,v);
divcircle = divF(x.^2 + y.^2 <= (u.^2 + v.^2));
flux_from_div = sum(divcircle)*h^2;