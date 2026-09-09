function f = ex2_fun(x,n)
    pn = zeros(size(x)); 
    for k = 0:n
        pn = pn + ((gamma(n+1))/((gamma(k+1)^2).*gamma(n-2*k+1))).*((2-x).^(n - 2*k)).*((x-1).^(k)).*(-1)^k;
    end
    f = sqrt(2*n + 1).*(x.^(-n-1)).*pn;
end