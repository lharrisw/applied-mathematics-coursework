%% Homework 2

% log-log plot

x = [1/10 1/100 1/1000 1/10000 1/100000 1/1000000 ...
    1/10000000 1/100000000 1/1000000000 1/10000000000];
y = [37 53 81 133 157 273 369 577 785 1173]; 

logx = log(x);
logy = log(y);
[p,S] = polyfit(logx,logy,1);
pp = polyval(p,logx);

figure
hold on
plot(logx,logy,'bo');
plot(logx,pp,'b-');
xlabel('$log(tol)$','interpreter','latex');
ylabel('$log(evals)$','interpreter','latex');
title(['Log of Function' ...
    ' Evaluations vs. Log of Tolerance'],'interpreter','latex');
legend('Data','Fit','interpreter','latex');
set(gca,'fontsize',20);
