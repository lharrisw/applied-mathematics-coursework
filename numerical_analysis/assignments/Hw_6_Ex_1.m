function Hw_6_Ex_1
    tic
    close all;
    clear;
    clc;
    
    N = 1000;
    p = initial_guess(N,1); % hopefully better guess
    x = optimize(p,eps('double'),20);
    u = energy(x)
    xx = linspace(-1,1,length(x));
    xx1 = linspace(-1,1,length(x));
%     poly(x);
%     txt1 = '$P_{1} = x$';
%     txt2 = '$P_{2} = x^2 - 1/3$';
%     txt3 = '$P_{3} = x^3 - 3/5x$';
%     txt4 = '$P_{4} = x^4 - 6/7x^2 + 3/35$';
%     txt5 = '$P_{5} = x^5 -10/9x^3 + 5/21$';
%     txt6 = '$P_{6} = x^6 - 15/11x^4 + 5/11x^2 - 5/231$';
%     txt = {txt1 txt2 txt3 txt4 txt5 txt6};
    
    % plots of the first 6 polynomials
%     figure
%     for k = 1:6
%         subplot(3,2,k)
%         hold on
%         plot(xx,LegendreP(k,xx),'b-');
%         plot(optimize(initial_guess(k,2),eps('double'),10),0,'ro');
%         title(txt{k},'interpreter','latex');
%     end

    % plots of the interlaced roots
%     figure
%     for j = 1:4
%         subplot(2,2,j)
%         hold on
%         plot(xx,LegendreP(j,xx),'b-',optimize(initial_guess(j,2),eps('double'),10),0,'bo');
%         plot(xx,LegendreP(j+1,xx),'r-',optimize(initial_guess(j+1,2),eps('double'),10),0,'ro');
%     end

   % distribution of the nodes
%    figure
%    hold on
%    my_hist(x);
    
   % 
% figure
% hold on
% plot(x,xx1,'b-');
% plot(x,asin(x)/pi+1/2,'r-');
% toc
end

function x = optimize(p,tol,MaxIter)
    opts.SYM = true;
    opts.POSDEF = true;
    x = p(:);
    o = ones(size(x));
    N = length(x);
    idx = (1:(N+1):N^2);
    c = 1/(2*pi*N^2);
    
    for k = 1:MaxIter
        
        % forming the gradient
        H = -1./(x*o' - o*x'); 
        y = 1./(1-x.^2);
        H(idx) = x.*y;
        
        g = sum(H,2); % gradient of u
        
        % forming the Hessian matrix
        H = -H.^2;
        H(idx) = 2*H(idx);
        H(idx) = y - sum(H,2);
        
        % solving linear system and 
        y = linsolve(c*H,c*g,opts);
        x = x - y;
        
        %fprintf('iter = %f; norm(y) = %f; cond(H) = %f\n',[k,norm(y,'inf'),cond(H)])
        
        % stopping point
        if norm(y,'inf') < tol
            break;
        end
    end
end

function x = initial_guess(N,MaxIter)
    x = 0; 
    for n = 1:N-1  
        x = .5*([-1;x] + [x;1]);
        x = optimize(x,eps('double'),MaxIter);
    end
end

function u = energy(p)
    x = p(:);
    N = length(x);
    c = 1/(2*pi*N^2); 
    o = ones(N,1);
    xx = tril(abs(x*o' - o*x'),-1);
    u = -c*(0.5*sum(log(1-x.^2)) + sum(log(xx(xx>0))));
end

function plt = my_hist(x)
    n = length(x);
    xx = linspace(-1,1,sqrt(n));
    dx = x(2) - x(1);
    yy = histc(x,xx);
    yy = yy/(numel(x)*dx);
    plt = bar(xx,yy);
end
