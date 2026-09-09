%% Exercise 1

n = 500;

x = randn(2,n);

plt1 = norm_plot(x,1);
plt2 = norm_plot(x,2);
plt3 = norm_plot(x,4);
plt4 = norm_plot(x,inf);
plt = {plt1 plt2 plt3 plt4};
txt1 = '1-Norm';
txt2 = '2-Norm';
txt3 = '3-Norm';
txt4 = 'Infinity-Norm';
txt = {txt1 txt2 txt3 txt4};

for k = 1:4
    subplot(2,2,k);
    plot(plt{k}(1,:),plt{k}(2,:),'b.'); 
    title(txt{k});
    xlim([-2 2]);
    ylim([-2 2]);
    axis equal
end

%% Exercise 7

x = linspace(-1,1,100);
xx = linspace(-1,1,8);

c = 4199*sqrt(1155)*[1 0 -36/19 0 378/323 0 -84/323 0 63/4199]/256;
root_s = roots(c)';

f = @(x) polyval(c,x);

y = @(x) 1./(1+25*x.^2);

p = polyfit(root_s,y(root_s),7);
pp = polyval(p,x);

p1 = polyfit(xx,y(xx),7);
pp1 = polyval(p1,x);


figure
hold on
plot(x,f(x));
plot(root_s,zeros(size(root_s)),'ro');
title('Plot of the Polynomial.')
set(gca,'fontsize',20);

figure
hold on
plot(x,y(x));
plot(x,pp,'k-');
plot(x,pp1,'r-');
legend('$\frac{1}{1+25x^{2}}$','Interpolant: Roots','Interpolants: Nodes','interpreter','latex');
set(gca,'fontsize',20);



%% Exercises 5 and 6

N = 4;

x = linspace(-1,1,1024);

% Minimizing 2-norm

r0 = -1 + 2*(1:N)/(N+1);
r = fminsearch(@two_norm,r0);
p = poly(r);

figure
plt = plot(x,abs(polyval(p,x)), 'b-',r,0,'ro',r0,0,'bo');

[p1,r1] = legendre_n(x,N);

% Minimizig 1-norm

rr0 = -1 + 2*(1:N)/(N+1);
rr = fminsearch(@one_norm,rr0);
pp = poly(rr);

figure
plt1 = plot(x,abs(polyval(pp,x)), 'b-',rr,0,'ro',rr0,0,'bo');

function f = one_norm(r)
    pp = poly(r);
    f = integral(@(x) abs(polyval(pp,x)),-1,1);
end

function f = two_norm(r)
    pp = poly(r);
    f = sqrt(integral(@(x)polyval(pp,x).^2,-1,1));
end