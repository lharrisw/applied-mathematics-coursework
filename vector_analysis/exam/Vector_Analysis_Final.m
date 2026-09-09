%% Exam 1 Exercise 1

% The Field
h = 0.01;
a = 1;
b = 1;
r = 1;
[x,y,z] = meshgrid(-5:h:5);
u = x;
v = y;
w = z;

% The faces of the cube
% Flux through face 1
u1 = x;
v1 = y;
w1 = 1;
t = 0:h:1
xt = t;
yt = t;
zt = t;
ut = interp3(x,y,z,u1,xt,yt,zt);
vt = interp3(x,y,z,v1,xt,yt,zt);
wt = interp3(x,y,z,w1,xt,yt,zt);

% Face 2
u2 = x;
v2 = 1;
w2 = z;

% Face 3
u3 = x;
v3 = y;
w3 = 0;

% Face 4
u4 = x;
v4 = 0;
w4 = y;

% Face 5
u5 = y;
v5 = z;
w5 = 1;

% Face 6
u6 = y;
v6 = z;
w6 = 0;

% Divergence through the interior
divF = divergence(x,y,u,v);
divcircle = divF((x-a).^2 + (y-b).^2 <= r^2);
flux_from_div = sum(divcircle)*h^2;