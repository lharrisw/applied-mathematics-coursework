%% Full Discretization of the Heat Equation: Euler's Method

N = 1000;
%[x,t] = meshgrid(linspace(0,5,N));

x = linspace(0,3,N);
t = linspace(0,1,N);

dx = length(x)/N;
dt = length(t)/N;
s = dt/dx^2;

u = zeros(length(x),length(t));

% Initial Condition
g = @humps;
u(:,1) = g(x);

% Boundary Condition
% u(1,:) = 0; % Dirichlet
% u(N,:) = 0; % Dirchlet

u(1,:) = 1; % Neumann
u(N,:) = 1; % Neumann

for tt = 1:N-1
    for xx = 2:N-2
        u(xx,tt+1) = u(xx,tt) + s*(u(xx+1,tt) - 2*u(xx,tt) + u(xx-1,tt));
    end
end


%% Wave Equation: Eulers Method

N = 1000;
%[x,t] = meshgrid(linspace(0,5,N));

x = linspace(0,3,N);
t = linspace(0,0.5,N);

dx = length(x)/N;
dt = length(t)/N;
s = (dt^2)/dx^2;

u = zeros(length(x),length(t));

% Boundary Condition
g = @humps;
u(:,1) = g(x);

% Initial Condition
% u(1,:) = 0; % Dirichlet
% u(N,:) = 0; % Dirchlet

for tt = 2:N-2
    for xx = 2:N-2
        u(xx,tt+1) = 2*u(xx,tt) - u(xx,tt-1) + s*(u(xx+1,tt) - 2*u(xx,tt) + u(xx-1,tt));
    end
end