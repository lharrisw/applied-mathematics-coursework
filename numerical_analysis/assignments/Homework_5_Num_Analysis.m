%% Exercise 1

omega = 100000;

f = @(s) omega*cos(omega*s)*0.5;
I = sin(omega);
tol = 1e-3;

Q_gauss = quad_gauss(f,100000,tol);
err_gauss = abs(I - Q_gauss);
Q_quadgk = quadgk(f,-1,1);
err_quadgk = abs(I - Q_quadgk);
Q_integral = integral(f,-1,1);
err_integral = abs(I - Q_integral);  