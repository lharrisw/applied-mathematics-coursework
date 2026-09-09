% Exercise 1

N = 100;
a = 2;
b = .1;
t = 1;
V = @(x,y) 10*(x<b & x>-b & y<b & y>-b) + 0*(x>b & x<-b & y>b & y<-b);
[x,y] = meshgrid(linspace(-a,a,N));

u = zeros(size(x));

% Sine Solution

% for m = 1:5
%     for n = 1:5
%         u = u + integral2(@(x,y) V(x,y).*sin(m*pi*x/(2*a)).*sin(n*pi*y/(2*a)),-a,a,-a,a)...
%             .*sin(m*pi*x/(2*a)).*sin(n*pi*y/(2*a)).*sin(sqrt(m^2 + n^2)*pi*t/(2*a))/...
%             integral2(@(x,y) (sin(m*pi*x/(2*a)).^2).*(sin(n*pi*y/(2*a)).^2),-a,a,-a,a)*(a/sqrt((pi*m)^2 + (pi*n)^2));
%     end
% end

% Cosine Solution

for m = 1:5
    for n = 1:5
        u = u + (a/(sqrt((m+1/2)^2 + (n+1/2)^2)*pi))*quad2d(@(xx,yy) V(xx,yy).*cos((m+1/2)*pi*xx/a).*cos((n+1/2)*pi*yy/a),-a,a,-a,a)...
            .*cos((m+1/2)*pi*x/a).*cos((n+1/2)*pi*y/a).*sin(sqrt((m+1/2)^2 + (n+1/2)^2)*pi*t/a)/...
            quad2d(@(xx,yy) (cos((m+1/2)*pi*xx/a).^2).*(cos((n+1/2)*pi*yy/a).^2),-a,a,-a,a);
    end
end

tt = linspace(0,10,1000);
du = zeros(size(tt));
x1 = 0.01;
y1 = 0.01;
ddu = zeros(size(tt));
xx1 = 0;
yy1 = 0;

% Sine Velocity: Always Zeros

% for n = 1:10
%     for m = 1:10
%         du = du + (sqrt(m^2 + n^2)*pi/a)*integral2(@(xx,yy) V(xx,yy).*sin(m*pi*xx/a).*sin(n*pi*yy/a),-a,a,-a,a)...
%             .*sin(m*pi*x1/a).*sin(n*pi*y1/a).*cos(sqrt(m^2 + n^2).*pi.*tt/a)/...
%             integral2(@(xx,yy) ((sin(m*pi*xx/a)).^2).*((sin(n*pi*yy/a)).^2),-a,a,-a,a);
%     end
% end

% Cosine Velocity

% for m = 1:15
%     for n = 1:15
%         ddu = ddu + sqrt((m+1/2)^2 + (n+1/2)^2)*quad2d(@(xx,yy) V(xx,yy).*cos((m+1/2)*pi*xx/a).*cos((n+1/2)*pi*yy/a),-a,a,-a,a)...
%             .*cos((m+1/2)*pi*xx1/a).*cos((n+1/2)*pi*yy1/a).*cos(sqrt((m+1/2)^2 + (n+1/2)^2)*pi*tt/a)/...
%             quad2d(@(xx,yy) (cos((m+1/2)*pi*xx/a).^2).*(cos((n+1/2)*pi*yy/a).^2),-a,a,-a,a);
%     end
% end


figure
mesh(x,y,u);
xlabel('$x$','interpreter','latex');
ylabel('$y$','interpreter','latex');
zlabel('$u(x,y,t)$','interpreter','latex');
title(sprintf('$u(x,y,1)$: $a$ = %1.0f, $b$ = %1.1f, $V$ = 10',[a b]),'interpreter','latex');
set(gca,'fontsize',18);
% 
% figure
% plot(tt,ddu);
% xlabel('Time (sec)','interpreter','latex');
% ylabel('$u_{t}(0,0,t)$','interpreter','latex');
% title(sprintf('Velocity: $a$ = %1.0f, $b$ = %1.1f, V = 10',[a b]),'interpreter','latex');
% set(gca,'fontsize',18);


%% Exercise 2

a = 2;
N = 100;
[x,y] = meshgrid(linspace(-a,a,N));

u1 = zeros(size(x));
u2 = zeros(size(x));
u3 = zeros(size(x));
u4 = zeros(size(x));

g1 = @(y) 1*(-a<y<a);
g2 = @(y) 1*(-a<y<a);
f1 = @(x) -1*(-a<x<a);
f2 = @(x) -1*(-a<x<a);

for n = 1:16    
    u1 = u1 + quadgk(@(yy) g1(yy).*sin(n*pi*yy/a),-a,a).*sinh(n*pi*(x-a)/a).*sin(n*pi*y/a)...
        /(a*sinh(-2*pi*n));
    u2 = u2 + quadgk(@(xx) f1(xx).*sin(n*pi*xx/a),-a,a).*sinh(n*pi*(y+a)/a).*sin(n*pi*x/a)...
        /(a*sinh(2*pi*n));
    u3 = u3 + quadgk(@(yy) g2(yy).*sin(n*pi*yy/a),-a,a).*sinh(n*pi*(x+a)/a).*sin(n*pi*y/a)...
        /(a*sinh(2*pi*n));
    u4 = u4 + quadgk(@(xx) f2(xx).*sin(n*pi*xx/a),-a,a).*sinh(n*pi*(y-a)/a).*sin(n*pi*x/a)...
        /(a*sinh(-2*pi*n));
end

uu = u1 + u2 + u3 + u4;

figure
mesh(x,y,uu);
xlabel('$x$','interpreter','latex');
ylabel('$y$','interpreter','latex');
zlabel('$u(x,y)$','interpreter','latex');
title('Solution to Laplace''s Equation','interpreter','latex');
set(gca,'fontsize',18);
axis square

%% Exercise 3

t = [0 0.1 0.5 1];

x = linspace(-5,5,100);

u = zeros(length(x),1);

for tt = 1:length(t)
    for j = 1:length(x)
        u(j,:) = integral2(@(w,y) exp(1i*w.*(x(j)-y))...
            .*exp(-(w.^2 + 1)*tt).*exp(-y.^2),-5,5,-5,5)/(2*pi);
    end
    hold on
    plot(x,u);
    xlabel('Rod Position','interpreter','latex');
    ylabel('$u(x,t)$','interpreter','latex');
    legend(sprintf('$t$ = %1.1f',t(1)),sprintf('$t$ = %1.1f',t(2))...
        ,sprintf('$t$ = %1.1f',t(3)),sprintf('$t$ = %1.1f',t(4)),'interpreter','latex');
    title('Temperature Over Time','interpreter','latex');
    set(gca,'fontsize',18)
end