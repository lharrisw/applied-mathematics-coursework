%% Exercise 3

% My Bessel 1

n = 100;
J = zeros(n,1);
z = linspace(0,10,n);

for k = 1:n
    j = 0;
    for l = 0:n
        j = j + (((-1)^l)*z(k)^(2*l+1))/(2^(2*l+1)*factorial(l)*factorial(l+1));
    end
    J(k,:) = j;
end

% Actual Besssel

y = besselj(1,z);

% Comparison

figure
hold on
plot(z,y,'ro');
plot(z,J,'b-');
xlabel('x','Interpreter',"latex");
ylabel('y','Interpreter',"latex");
title('Bessel Function of the First Kind, Order One');
legend('Actual Bessel Function','My Bessel Function');
set(gca,'fontsize',20);

%% Exercise 4

% Orthogonality of Bessel Functions

n = 100;
z = linspace(0,10,100);
J0 = @(z) besselj(0,z);

% Generating a list of zeros

for k = 1:20
    zn(k) = fzero(J0,k);
end

% Obtaining unique zeros

zn = uniquetol(zn,eps(1.0)); 
zn = zn(1:4);
prod1 = zeros(4,4);

% Checking the orthogonality of the bessel function

for ii = 1:length(zn)
    for jj = 1:length(zn)
        prod1(ii,jj) = quadgk(@(r)besselj(0,zn(ii)*r).*besselj(0,zn(jj)*r).*r,0,1);
    end
end

prod2 = zeros(length(zn),1);

% Integrating the first four terms: n = 0,1,2,3

for kk = 1:length(zn)
    prod2(kk,:) = quadgk(@(r)besselj(0,zn(kk)*r).^(2).*r,0,1);
end