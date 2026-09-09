function Q = quad_gauss_exam(f,n,tol)

x0 = cos(pi*(0.75 + (0:n-1))/(n + 0.5));
y0 = cos(pi*(0.75 + (0:n-1))/(n + 0.5));

for k = 1:1000
    [p,dp] = legendrep(n,x0);
    x0 = x0 - (p./dp).*(x0.^2 - 1);
    if norm(p,inf) < tol
        break
    end
end

for k = 1:1000
    [p,dp] = legendrep(n,y0);
    y0 = y0 - (p./dp).*(y0.^2 - 1);
    if norm(p,inf) < tol
        break
    end
end

w1 = 2./(dp.^2).*(1 - x0.^2);
w2 = 2./(dp.^2).*(1 - y0.^2);

o = ones(size(x0));
p = (x0'*o).^2 + (o'*y0).^2;

Q = w1*f(p)*w2';
end

function [p,dp] = legendrep(n,x)

if n == 0
    p = ones(size(x));
    dp = zeros(size(x));
elseif n == 1
    p = x;
    dp = ones(size(x));
else
    p1 = x;
    p2 = ones(size(x));
    for m = 2:n
        p = ((2*m - 1)*(x.*p1) - (m-1)*p2)/m;
        p2 = p1;
        p1 = p;
    end
    dp = n.*(x.*p - p2);
end    
end