function Hw_6_Ex_3
    close all;
    clc;
    
    N = 2;
%     x = linspace(-1,1,N);
%     xx = linspace(1,1,N);
    g = initial_guess(N,eps('double'),1); % much better guess
    t = optimize(g,eps('double'),25);
%     u = energy(t);
%     p = 6*poly(t);
%     pl = [cos(t) sin(t)]';
    
    % plot of first seven polynomials
%     figure
%     subplot(4,2,[1 2])
%     r = optimize(initial_guess(0,eps('double'),2),eps("double"),5);
%     rr = linspace(0,pi);
%     p = poly(r);
%     hold on
%     plot(optimize(initial_guess(0,eps('double'),2),eps('double'),5),0,'ro');
%     plot(rr,polyval(p,rr),'b-');
%     xlabel('$x$','interpreter','latex');
%     ylabel('$y$','interpreter','latex');
%     set(gca,'fontsize',15);
%     for k = 0:5
%         r = optimize(initial_guess(k+2,eps('double'),2),eps("double"),5);
%         rr = linspace(0,pi);
%         p = poly(r);
%         subplot(4,2,k+3)
%         hold on
%         plot(optimize(initial_guess(k+2,eps('double'),2),eps('double'),5),0,'ro');
%         plot(rr,polyval(p,rr),'b-');
%         xlabel('$x$','interpreter','latex');
%         ylabel('$y$','interpreter','latex');
%         set(gca,'fontsize',15);
%     end
    
    % plot of interlaced nodes
%     figure
%     for j = 1:4
%         r = optimize(initial_guess(j,eps("double"),2),eps("double"),2);
%         r1 = optimize(initial_guess(j+1,eps("double"),2),eps("double"),2);
%         rr = linspace(0,pi);
%         p = poly(r);
%         p1 = poly(r1);
%         subplot(2,2,j)
%         hold on
%         plot(rr,polyval(p,rr),'b-',optimize(initial_guess(j,eps("double"),2),eps('double'),2),0,'bo');
%         plot(rr,polyval(p1,rr),'r-',optimize(initial_guess(j+1,eps("double"),2),eps('double'),2),0,'ro');
%         xlabel('$x$','interpreter','latex');
%         ylabel('$y$','interpreter','latex');
%         set(gca,'fontsize',20);
%     end
    
    % plot of the wire with the nodes
%     figure
%     hold on
%     plot(pl(1,:),pl(2,:),'ro');
%     plot(pl(1,:),pl(2,:),'b-')
%     xlabel('$x$','interpreter','latex');
%     ylabel('$y$','interpreter','latex');
%     set(gca,'fontsize',20);
%     title('Plot of the $p_{i}$ on the Wire','interpreter','latex');
%     axis equal
    
    % plot of the distribuion of the nodes
%     figure
%     my_hist(t);
%     xlabel('$x$','interpreter','latex');
%     ylabel('$y$','interpreter','latex');
%     title('Distribution of Nodes','interpreter','latex');
%     set(gca,'fontsize',20);
%     
%     figure
%     hold on
%     plot((t)/max(t),x,'b-');
%     plot(x/2+0.5,2*asin(x)/pi,'r-');
%     xlabel('$x$','interpreter','latex');
%     ylabel('$y$','interpreter','latex');
%     title('The CDF of the Nodes','interpreter','latex');
%     legend('CDF of the Nodes','Shifted \& Scaled Inverse Sine Law','interpreter','latex');
%     set(gca,'fontsize',20);
 
end

function t = optimize(p,tol,MaxIter)
    opts.SYM = true;
    opts.POSDEF = true;
    t = p(:);
    o = ones(size(t));
    N = length(t);
    idx = (1:(N+1):N^2);
    
    for k = 1:MaxIter
        
        % forming the gradient
        tt = 0.5*(t*o'-o*t');
        H = 0.5*cot(tt);
        H(idx) = cot(t);
        g = -sum(H,2);
        
        % forming the Hessian
        H = H.^2;
        H(idx) = 0;
        HH = 0.25 + H;
        HH(idx) = -csc(t).^2 - sum(H,2);
        HH = -HH;
        
        % solving linear system and 
        y = linsolve(HH,g,opts);
        t = t - y;
        
        %fprintf('iter = %f; norm(y) = %f; cond(H) = %f\n',[k,norm(y,'inf'),cond(H)])
        
        % stopping point
        if norm(y,'inf') < tol
            break;
        end 
    end
end

function t = initial_guess(N,tol,MaxIter)
    t = pi/2;    
    for n = 1:N-1
        t = 0.5*([0;t] + [t;pi]);
        t = optimize(t,tol,MaxIter); 
    end
end

function u = energy(t)
    % output of method 1 = output of method 2
    
    % method 1
    N = length(t);
    s1 = log(sin(t));
    tt = t(:);
    c = 1/(2*pi*N^2);
    o = ones(N,1);
    s2 = 0.5*abs(tt*o' - o*tt');
    s2 = tril(s2);
    s2 = log(sin(s2(s2>0)));
    u = -c*(2*N*log(2) + sum(s1) + sum(s2));

    % method 2
%     p0 = [cos(0);sin(0)];
%     pN1 = [cos(pi);sin(pi)];
%     p = [cos(t);sin(t)];
%     s1 = log(vecnorm(p0 - p));
%     s2 = log(vecnorm(pN1 - p));
%     tt = t(:);
%     N = length(t);
%     o = ones(N,1);
%     c = 1/(2*pi*N^2);
%     s3 = 0.5*abs(tt*o' - o*tt');
%     s3 = tril(s3);
%     s3 = (log(sin(s3(s3>0))));    
%     u = -c*(N*log(2) + sum(s1) + sum(s2) + sum(s3));
end

function plt = my_hist(x)
    n = length(x);
    xx = linspace(0,pi,sqrt(n));
    dx = x(2) - x(1);
    yy = histc(x,xx);
    yy = yy/(numel(x)*dx);
    plt = bar(xx,yy);
end