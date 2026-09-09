function Hw_6_Ex_2
tic
close all;
clear;
clc;

N = 500;
p = initial_guess(N,2); % hopefully better guess
x = optimize(p,eps('double'),20);
u = energy(x);
%xx = linspace(-1,1,N);
    
% plots of the first 7 polynomials
figure
subplot(4,2,[1 2])
r = optimize(initial_guess(0,2),eps("double"),20);
rr = linspace(-1,1);
p = poly(r);
hold on
plot(optimize(initial_guess(0,2),eps('double'),10),0,'ro');
plot(rr,polyval(p,rr),'b-');
xlabel('$x$','interpreter','latex');
ylabel('$y$','interpreter','latex');
set(gca,'fontsize',15);
for k = 0:5
    r = optimize(initial_guess(k+2,2),eps("double"),20);
    rr = linspace(-1,1);
    p = poly(r);
    subplot(4,2,k+3)
    hold on
    plot(optimize(initial_guess(k+2,2),eps('double'),20),0,'ro');
    plot(rr,polyval(p,rr),'b-');
    xlabel('$x$','interpreter','latex');
    ylabel('$y$','interpreter','latex');
    ylim([-0.5 0.5]);
    set(gca,'fontsize',15);
end

%     plots of the interlaced roots
    figure
    title('Interlaced Nodes');
    for j = 1:4
        r = optimize(initial_guess(j,2),eps("double"),10);
        r1 = optimize(initial_guess(j+1,2),eps("double"),10);
        rr = linspace(-1,1);
        p = poly(r);
        p1 = poly(r1);
        subplot(2,2,j)
        hold on
        plot(rr,polyval(p,rr),'b-',optimize(initial_guess(j,2),eps('double'),10),0,'bo');
        plot(rr,polyval(p1,rr),'r-',optimize(initial_guess(j+1,2),eps('double'),10),0,'ro');
        xlabel('$x$','interpreter','latex');
        ylabel('$y$','interpreter','latex');
        set(gca,'fontsize',20);
    end

%     distribution of the nodes
    figure
    hold on
    my_hist(x);
    xlabel('$x$','interpreter','latex');
    ylabel('$y$','interpreter','latex');
    title('Distribution of Nodes','interpreter','latex');
    set(gca,'fontsize',20);
     
%     plot of CDF 
    figure
    hold on
    plot(x,xx,'b-');
    xlabel('$x$','interpreter','latex');
    ylabel('$y$','interpreter','latex');
    title('The CDF of the Nodes','interpreter','latex');
    set(gca,'fontsize',20);
     toc
end

function x = optimize(p,tol,MaxIter)
    opts.SYM = true;
    opts.POSDEF = true;
    x = p(:);
    o = ones(size(x));
    N = length(x);
    idx = (1:(N+1):N^2);
    
    for k = 1:MaxIter
        
        % forming the gradient
        xx = x*o'-o*x';
        H = xx.^-2;
        H = -tril(H,-1) + triu(H,1);
        y = 4./(1-x.^2).^2;
        H(idx) = x.*y;
        
        g = sum(H,2);
        
        % Forming the Hessian
        H = 2*H./xx;
        H(idx) = -y.^2.*x.^2.*(1-x.^2);
        H(idx) = y - sum(H,2);
        
        % solving linear system and 
        y = linsolve(H,g,opts);
        x = x - y;
        
        %fprintf('iter = %f; norm(y) = %f; cond(H) = %f\n',[k,norm(y,'inf'),cond(H)])
        
        % stopping point
        if norm(y,'inf') < tol
            break;
        end 
    end
end

function x = initial_guess(N,MaxIter)
    opts.SYM = true;
    opts.POSDEF = true;
    x = 0;
    for k = 1:N-1
        x = 0.5*([-1;x] + [x;1]);
        x = optimize(x,eps('double'),MaxIter);
%         o = ones(k+1,1);
%         idx = (1:(k+2):(k+1)^2);
%         xx = x*o'-o*x';
%         H = xx.^-2;
%         H = -tril(H,-1) + triu(H,1);
%         y = 4./(1-x.^2).^2;
%         H(idx) = y.*x;
%         g = sum(H,2);
%         H = 2*H./xx;
%         H(idx) = -y.^2.*x.^2.*(1-x.^2);
%         H(idx) = y - sum(H,2);
%         y = linsolve(H,g,opts);
%         x = x - y;
    end
end

function u = energy(p)
    x = p(:);
    N = length(x);
    c = 1/(4*pi*N^2); 
    o = ones(N,1);
    xx = 1./abs(x*o' - o*x');
    xx = tril(xx,-1);
    u = c*(sum(2./(1-x.^2)) + sum(xx(:)));
end

% Justin Energy Function 
% function u = energy1(x0)
%     x = x0(:);
%     N = length(x);
%     ind = (1:(N+1):N^2)';
%     o = ones(size(x));
%     xx = x*o'-o*x';
%     A = 1./abs(xx);
%     A(ind) = 4./(1-x.^2);
%     u = sum(A(:))/(8*pi*N^2);
% end

function plt = my_hist(x)
    n = length(x);
    xx = linspace(-1,1,sqrt(n));
    dx = x(2) - x(1);
    yy = histc(x,xx);
    yy = yy/(numel(x)*dx);
    plt = bar(xx,yy);
end