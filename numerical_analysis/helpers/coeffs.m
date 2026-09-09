function [c,gram] = coeffs(N)
    gram = zeros(N,N);
    for ii = 1:N
        for jj = 1:N
            gram(ii,jj) = integral(@(x) LegendreP(ii-1,x).*LegendreP(jj-1,x),-1,1);
        end
    end 
    b = zeros(1,N);
    for k = 1:N
        b(k) = integral(@(x) exp(x).*LegendreP(k-1,x),-1,1);
    end
    c = pinv(gram)*b';
    %c = gram\b';
end