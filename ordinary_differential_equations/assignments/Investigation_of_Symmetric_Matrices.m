%% Investigation of Symmetric Matrices

n = 100;
m = 100;

% Lower Triangular Matrix with Ones
% for j = 1:n;
%     for k = 1:j;
%         A(j,k) = 1;
%     end
% end
% 
% for l = 1:n;
%     for m = 1:l;
%         B(l,m) = n-l+m;
%         B(m,l) = B(l,m);
%     end
%     B(l,l) = n;
% end

A = toeplitz([1 -2 1 100 0 0]);
[V,D] = eig(A);

B = randn(m,n);
[v,d] = eig(B);


figure
hold on
plot(diag(D),'b.');

figure 
hold on
plot(diag(d),'k.');

%% Matrix vector ODE

t = linspace(0,10,1000);
A = toeplitz([-2 1]);
dxdt = @(t,x) [x(2); -4*x(2) - 3*x(1)];
x0 = [1;0];
[t,x] = ode45(dxdt,[0 10],x0);
x1 = -.5*exp(-3.*t) + 1.5*exp(-t);

figure
hold on
plot(t,x(:,1),'b-');
plot(t,x1,'ro');

%% solving the Heat equation

N = 1000;
L = 10;
h = L/(N+1);
x = (1:N)*h;
A = toeplitz([-2 1 zeros(1,N-2)]);
A(1,1) = -1;
A(N,N) = -1;

[V,D] = eig(A);
lambda = diag(D)/(h^2);
f = @humps;
u0 = f(x)';
t = 0.0001;
C = V\u0;
u = V*(exp(lambda*t).*C);
tt = linspace(0,1,N);

figure
hold on
p = plot(x,u0,'b-');
plot(x,u,'r-');


for t = tt
    u = V*(exp(lambda*t).*C);
    if ishandle(p)
        set(p,'ydata',u);
        title(sprintf('t=%1.4f',t),'FontSize',20);
        pause(0.0005);
    end
end

%% TST Eigenvalue & Eigenvector derivation

N = 100;
L = 1;

h = L/(N+1);
x = (1:N)*h;

A = toeplitz([2 -1 zeros(1,N-2)]);
[V,D] = eig(A);

[lambda,l] = sort(abs(diag(D)));
V = V(:,l);

% for k = 1:9
%     subplot(3,3,k);
%     hold on
%     plot(x,sign(V(1,k)).*V(:,k)/abs(max(V(:,k))),'r.');
%     plot(x,sin(k*x));
% end

for k = 1:9
    subplot(3,3,k);
    hold on
    plot(x,(1/(pi*k))*asin(sign(V(1,k)).*V(:,k)*sqrt((N+1)/2)),'r.');
    %plot(x,sqrt(2/(N + 1))*sin(x*k*pi));
    plot(x,x);
end