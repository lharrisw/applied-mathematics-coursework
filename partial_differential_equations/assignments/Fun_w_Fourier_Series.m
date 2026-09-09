%% Fun with Fourier Series

N = 250;
x = linspace(0,2*pi,N);

%y = @(x) exp(x);
%y = @sqrt;
%y = @square;
%y = @(x) 4*x + 3*log(x);
%y = @sawtooth;
y = @sinc;

yy = zeros(size(x));

% Fourier Series of the function y

for n = 1:5
    yy = yy + (1/integral(@(t) sin(n*t).^2,0,2*pi))*integral(@(t) y(t).*sin(n*t),0,2*pi).*sin(n*x)...
        + (1/integral(@(t) cos(n*t).^2,0,2*pi))*integral(@(t) y(t).*cos(n*t),0,2*pi).*cos(n*x);
end

yy = yy + integral(@(t) y(t),0,2*pi)/(2*pi); % adding zeroth term in fourier sum

figure
hold on
plot(x,y(x),'r-');
plot(x,yy,'b-');