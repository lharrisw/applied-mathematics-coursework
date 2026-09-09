%% Exercise 1

w = linspace(0,2);
a = 2.5;
b = 1.1;
c = .75;
G = 1./sqrt((c - a*(w.^2)).^2 + (b*w).^2);
dG = diff(G);
dw = diff(w);
dGdw = dG./dw;
[yvalue,xvalue] = max(G);
dGdw1 = -((2*(a^2)*(w.^3) + (b^2)*(w) - 2*a*c*w)./(((a^2)*(w.^4) + (b^2)*(w.^2) - 2*a*c*(w.^2) + (c^2)).^(3/2)));
wres = (1/a)*sqrt((2*a*c - b^2)/2);

figure
hold on
plot(w,G,'b-');
plot(w(1:end-1),dGdw,'mo');
plot(w,dGdw1,'r-');
plot(wres,yvalue,'ko');
xlabel('\omega');
ylabel('G(\omega)');
legend('G(\omega)', 'Numerical Derivative of G', 'Symbolic Derivative of G','(\omega_{res}, G(\omega_{res}))');
%% Exercise 2

R = .1;
L = 1;
C = .1;
a = .33333;
b = .66666;
w = linspace(0,10);

%% Exercise 3

t = linspace(0,2*pi);
n = 1:2:(2*size(t'));
R = 0.1;
L = 1;
C = 10;
y = asin(sin(t));
b_n = 2*sin(n*pi*.5)./(pi*n.^2) - 2*sin(3*n*pi*.5)./(pi*n.^2);
f = b_n*sin(n'*t);
%p1 = -L*n.^2 + R*i*n + (1/C);
%p2 = -L*n.^2 - R*i*n + (1/C);
%xp = (b_n./(2*i))'*(exp(i*n'*t)/p1 - exp(-i*n'*t)/p2);
dxdt = @(t,x) [x(2); (f - R*x(2) - (1/C)*x(1)/L)];
x0 = [0;0];
T = 2*pi;
[t,x] = ode45(dxdt,[0 T],x0);

figure
hold on
plot(t,x(:,1),'b-');
%plot(t,y,'r-');
%plot(t,f,'b-');