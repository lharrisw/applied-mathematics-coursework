%% PDE Exam Maple solution

N = 50;
% a = 2;
% a = 4;
% a = 8;
% a = 16;

b = 1;

[r,theta] = meshgrid(linspace(0,a,N),linspace(0,2*pi,N));

c2 = 1;
C2 = 1;
C3 = 3;
C4 = randn;

% u = @(r,theta) (c2*log(r).^2)/2 + log(b)*c2*log(r) ...
%     + C2 + 0.*(r < b) + (r >= b).*((r.^2)/4 - (b^2)/4 + (log(b)*b^2)/2 - (log(r)*b^2)/2)...
%     + (c2*theta.^2)/2 + C3*theta + C4;

u = @(r,theta)C2 + (c2*log(b) - (b^2)/2).*log(r) - (c2*log(r).^2)/2 ...
    + 0.25*(r.^2).*(r < b) +  (-b^2*(-1 + 2*log(b) - 2*log(r))).*(r >= b) ...
    + (c2*theta.^2)/2 + C3*theta + C4;

figure
mesh(r.*cos(theta),r.*sin(theta),u(r,theta)); 
xlabel('$x$','Interpreter','latex');
ylabel('$y$','Interpreter','latex');
zlabel('$u(r,\theta)$','Interpreter','latex');
title(sprintf('Poisson Equation Solution: $a$ = %1.0f',a),'interpreter','latex');
set(gca,'fontsize',18);