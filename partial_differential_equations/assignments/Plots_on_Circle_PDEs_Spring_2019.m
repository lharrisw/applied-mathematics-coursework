%% Continuos functions
clear
close all
clc

t = linspace(-pi,pi,1000);
x = cos(t);
y = sin(t);
z = (t.^2 - cos(t));
dz = diff(z);

plot3(x,y,z);
plot3(x(1:end-1),y(1:end-1),dz);
%plot3(x(1:end-2),y(1:end-2),ddz);