function [lambda,lambda1] = dir_eig(N,n)
lambda = 2*sin((pi*(N - 2*n + 1))/(2*N + 2)) - 2; % The derived formula
lambda1 = 2*cos((pi*n)/(N + 1)) - 2; % The actual formula
end