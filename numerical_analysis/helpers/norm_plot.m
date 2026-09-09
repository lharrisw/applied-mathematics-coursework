function n = norm_plot(x,p)
n = zeros(size(x));
for c = 1:length(x)
    n(:,c) = x(:,c)/norm(x(:,c),p); 
end
end