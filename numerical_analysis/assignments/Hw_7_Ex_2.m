function Hw_7_Ex_2
    clear;
    close all;
    clc;
    
    f = @(t,s) 0.*t.*s;
    lambda = -1;
    f1 = @(t,s) lambda*s;
    h = [0.1 0.05 0.025];
    
    t1 = 0:h(3):20;
    
    figure
    plot(t1,exp(lambda*t1),'b-');
    for k = 1:3
        hold on
        [y1,t1] = ode_solver(f1,0,20,h(k),1,exp(lambda*h(k)));
        plot(t1,y1,'-.');
        ylim([-2.5 2.5]);
        xlim([0 2.5]);
    end
    
    [y,t] = ode_solver(f,0,20,h(2),1,1);
    
    figure
    hold on
    plot(t,ones(size(t)),'b-')
    plot(t,y,'r-.')
    
end

function [y,t] = ode_solver(rhs,t0,tf,h,y0,y1)
    t = t0:h:tf;
    y = zeros(size(t));
    y(1) = y0;
    y(2) = y1;
    
    for k = 1:length(t)-2
        ynew = y(k+1) + rhs(t(k+1),y(k+1))*h;
        for j = 1:100
            ynew = 3*y(k+1) - 2*y(k)...
            + h*(13*rhs(t(k+2),ynew)/12 - 5*rhs(t(k+1),y(k+1))/3 - 5*rhs(t(k),y(k))/12);
        end
        y(k+2) = 3*y(k+1) - 2*y(k)...
            + h*(13*rhs(t(k+2),ynew)/12 - 5*rhs(t(k+1),y(k+1))/3 - 5*rhs(t(k),y(k))/12);
    end
end