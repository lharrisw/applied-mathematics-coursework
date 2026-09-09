%% Exercise 1
clear
close all
clc

n = 1000;
A = ones(1,n);

for i = 1:n
    a = 0;
    for j = 1:i
        u = randn(2,1);
        v = randn(2,1);
        a = a + sqrt(((norm(u)^2)*(norm(v)^2) - (dot(u,v)^2)));
    end
    A(i) = a/i;
end

figure
hold on
plot(1:n,A);
title('R^{2}');
xlabel('Number of Parallelograms');
ylabel('Average Area');


%% Test 1 Exercise 1
clear
close all
clc

n = 1000;
A = ones(1,n);

for i = 1:n
    a = 0;
    for j = 1:i
        u = randn(3,1);
        v = randn(3,1);
        a = a + sqrt(((norm(u)^2)*(norm(v)^2) - (dot(u,v)^2)));
    end
    A(i) = a/i;
end

figure
hold on
plot(1:n,A);
title('R^{3}')
xlabel('Number of Parallelograms');
ylabel('Average Area');
ylim([0 5]);

%% Test 1 Exercise 1
clear
close all
clc

n = 1000;
A = ones(1,n);

for i = 1:n
    a = 0;
    for j = 1:i
        u = randn(4,1);
        v = randn(4,1);
        a = a + sqrt(((norm(u)^2)*(norm(v)^2) - (dot(u,v)^2)));
    end
    A(i) = a/i;
end

figure
hold on
plot(1:n,A);
title('R^{4}');
xlabel('Number of Parallelograms');
ylabel('Average Area');
ylim([0 5]);

%% Test 3 Exercise 1
clear
close all
clc

n = 1000;
A = ones(1,n);

for i = 1:n
    a = 0;
    for j = 1:i
        u = randn(5,1);
        v = randn(5,1);
        a = a + sqrt(((norm(u)^2)*(norm(v)^2) - (dot(u,v)^2)));
    end
    A(i) = a/i;
end

figure
hold on
plot(1:n,A);
title('R^{5}');
xlabel('Number of Parallelograms');
ylabel('Average Area');
ylim([0 6]);

%% Test 4 Exercise 1
clear
close all
clc

n = 1000;
A = ones(1,n);

for i = 1:n
    a = 0;
    for j = 1:i
        u = randn(6,1);
        v = randn(6,1);
        a = a + sqrt(((norm(u)^2)*(norm(v)^2) - (dot(u,v)^2)));
    end
    A(i) = a/i;
end

figure
hold on
plot(1:n,A);
title('R^{6}')
xlabel('Number of Parallelograms');
ylabel('Average Area');
ylim([0 7]);

%% Test 5 Exercise 1
clear
close all
clc

n = 1000;
A = ones(1,n);

for i = 1:n
    a = 0;
    for j = 1:i
        u = randn(7,1);
        v = randn(7,1);
        a = a + sqrt(((norm(u)^2)*(norm(v)^2) - (dot(u,v)^2)));
    end
    A(i) = a/i;
end

figure
hold on
plot(1:n,A);
xlabel('Number of Parallelograms');
ylabel('Average Area');
title('R^{7}');
ylim([0 8]);

%% Exercise 3
clear
close all
clc

[x,y] = meshgrid(-1:0.1:1,-1:0.1:1);
w = 2;
u = -y;
v = x;

figure
quiver(x,y,u,v);
axis equal

%% Exercise 4
clear
close all
clc

h = [0 0 1];
u = randn(1,3);
u = u./norm(u);

v = cross(h,u);
v = v./norm(v);

n = cross(u,v);
n = n./norm(n);

w = 0.2;
r = 1;
[x,y,z] = meshgrid(-r-1:r/2:r+1);

f1 = -w.*r.*y;

A = [v' n' u'];

for i = 1:length(f1);
    for j = 1:length(f1);
        for k = 1:length(f1);
            f = (A\[x(i,j,k);y(i,j,k);z(i,j,k)])';
            dv(i,j,k) = f(1,1);
            dn(i,j,k) = f(1,2);
            du(i,j,k) = f(1,3);
        end
    end
end

f1 = -w.*r.*dn;
f2 = w.*r.*dv;
f3 = 0*f1;

for i = 1:length(f1)
    for j = 1:length(f1)
        for k = 1:length(f1)
            f = (A*[f1(i,j,k);f2(i,j,k);f3(i,j,k)])';
            fv(i,j,k) = f(1,1);
            fn(i,j,k) = f(1,2);
            fu(i,j,k) = f(1,3);
        end
    end
end

[a,b,c] = sphere;
surf(a*r,b*r,c*r);
hold on
quiver3(x,y,z,fv,fn,fu);
plot3([0 u(1)*(r+1)],[0 u(2)*(r+1)],[0 u(3)*(r+1)]);
axis equal

%% Exercise 5
clear
close all
clc

% The Field
h = 0.01;
a = 1;
b = 1;
r = 1;
[x,y] = meshgrid(-5:h:5);
z = peaks(x,y);
[u,v] = gradient(z,h,h);

% Flux across boundary
t = 0:h:2*pi;
xt = a + r*cos(t);
yt = b + r*sin(t);
ut = interp2(x,y,u,xt,yt);
vt = interp2(x,y,v,xt,yt);
Fnorm = ut.*cos(t) + vt.*sin(t);
flux = sum(Fnorm*h*r);
disp(flux)

% Divergence through the interior
divF = divergence(x,y,u,v);
divcircle = divF((x-a).^2 + (y-b).^2 <= r^2);
flux_from_div = sum(divcircle)*h^2;
disp(flux_from_div)

%% Exam 1 Exercise 1c
clear
close all
clc

% The Field
h = 0.01;
a = 1;
b = 1;
r = 1;
[x,y] = meshgrid(-5:h:5);
u = x.^3;
v = y.^3;

% Flux across boundary
t = 0:h:2*pi;
xt = a + r*cos(t);
yt = b + r*sin(t);
ut = interp2(x,y,u,xt,yt);
vt = interp2(x,y,v,xt,yt);
Fnorm = ut.*cos(t) + vt.*sin(t);
flux = sum(Fnorm*h*r);

% Divergence through the interior
divF = divergence(x,y,u,v);
divcircle = divF((x-a).^2 + (y-b).^2 <= r^2);
flux_from_div = sum(divcircle)*h^2;


%% Exercise 6
clear
close all
clc

x0 = -1;
y0 = -1;
r = 1;
x = @(t) x0 + r.*cos(t);
y = @(t) y0 + r.*sin(t);
n = @(t) [cos(t);sin(t)];
F = @(t) [x(t)./(x(t).^2 + y(t).^2); y(t)./(x(t).^2 + y(t).^2)];
Fn = @(t) F(t)'*n(t)*r;

area = integral(Fn,0,2*pi,'ArrayValued',true);

[a,b] = meshgrid(-3:.25:3,-3:.25:3);
u = a./(a.^2 + b.^2);
v = b./(a.^2 + b.^2);

figure
hold on
quiver(a,b,u,v);
plot(x(0:0.01:2*pi),y(0:0.01:2*pi));
axis equal       