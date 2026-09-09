function [component] = dir_vec(N,n,k)
f = (2*k*n)/(N - 2*n + 1);
component = sqrt(2/(N+1))*sin((pi*(N - 2*n + 1)*f)/(2*N + 2));
end