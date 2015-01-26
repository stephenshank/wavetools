clear all
close all
clc

f = 10;			% frequency
n = 30*f;		% number of grid points
h = 1/(n+1);
omega = 2*pi*f;	% angular frequency
g = @(x) 1i*exp(1i*omega*abs(x-.5))/(2*omega);

x = .5;			% coordinate of point source
m = ones(n+2,1);
A = helmholtz_fos_1d(m,f);
b = zeros(n+2,1);
b(round(x/h)) = 1/h;
u = A\b;
x=linspace(0,1,n+2);
subplot(1,2,1)
plot(x,real(u),'-rx',x,real(g(x)),'-bo')
subplot(1,2,2)
plot(x,imag(u),'-rx',x,imag(g(x)),'-bo')