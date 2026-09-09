function y = InvSqrt(x)
    tol = 1e-2;
    y = newt_boi(x,tol);
    y = abs(y);
end

% function xx = newt_boi(x,tol)
%     xx = 1.5 - 0.5*x; 
%     while any(xx.^-2 - x > tol)
%         xx = 1.5*xx - 0.5*x.*xx.*3;
%     end
% end

function xx = newt_boi(x,tol)
    xx = sqrt(x); 
    while any(xx.^2 - 1./x > tol)
        xx = 0.5*xx + 0.5*(1./x).*(1./xx);
    end
end