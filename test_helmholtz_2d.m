% Demonstrate usage of the Helmholtz and domain classes on common community
% models.
%% Clear and close everything
clear all
close all
clc

%% Key parameters for setup
f = 20;			% frequency
n = 10*f;		% number of grid points
omega = 2*pi*f;	% angular frequency
x = .2;			% x coordinate of point source
y = .8;			% x coordinate of point source
pml = struct('width',.1,'intensity',5e4);

%% Free space solution
figure
dom = domain([0 1 0 1],[n n]);
m = ones(dom.N,1);
A = helmholtz_2d(m,f,dom,pml);
b = dom.pt_src(x,y);
subplot(2,4,5)
dom.imagesc(real(A\b))
title('Discrete approximation')
cvec = caxis;

subplot(2,4,1)
Utrue = flipud(1i*besselh(0,omega*abs(dom.X+1i*dom.Y-(x+1i*y)))/4);
imagesc(real(Utrue));
title('Free space solution')
caxis(cvec)

%% Shepp-Logan phantom
dom = domain([0 1 0 1],[n n]);
C = phantom2(dom);
m = dom.mat2vec(1./C.^2);
A = helmholtz_2d(m,f,dom,pml);
b = dom.pt_src(x,y);
subplot(2,4,6)
dom.imagesc(real(A\b))
title('Discrete approximation')
caxis(cvec)

subplot(2,4,2)
dom.imagesc(C);
hold on
dom.plot(x,y,'r.','MarkerSize',10);
hold off
title('Shepp-Logan phantom')
caxis([1 2])

%% Marmousi model
subplot(2,4,7:8)
dom = domain([0 384/122 0 1],[ceil(384/122*n) n]);
C = marmousi(dom);
m = dom.mat2vec(1./C.^2);
A = helmholtz_2d(m,f,dom,pml);
b = dom.pt_src(x,y);
dom.imagesc(real(A\b))
title('Discrete approximation')
caxis(cvec)

subplot(2,4,3:4)
dom.imagesc(C);
hold on
dom.plot(x,y,'r.','MarkerSize',10);
hold off
title('Marmousi')
caxis([1 2])