clear all
close all
clc

maxit=200;
freq.min=1;
freq.max=5;
freq.n=6;

receiver.left=[.05 .15];
source.left=[.2 .3];
source.right=[.7 .8];
receiver.right=[.85 .95];

source.n=6;
receiver.n=6;

c=@(x) 1+.2*exp(-100*(x-.5).^2);
n=20*freq.max;
h=1/(n+1);
x=h:h:1-h;
noise=1e-12;
method='lbfgs';

sl=source.left;
rl=receiver.left;
sr=source.right;
rr=receiver.right;
ns=source.n;
nr=receiver.n;
xs_ind=[round(linspace(sl(1),sl(2),ns/2)/h) round(linspace(sr(1),sr(2),ns/2)/h)];
xr_ind=[round(linspace(rl(1),rl(2),nr/2)/h) round(linspace(rr(1),rr(2),nr/2)/h)];

figure
plot(x(xs_ind),ones(size(xs_ind)),'rx',x(xr_ind),ones(size(xr_ind)),'bo',x,c(x),'-k')
legend('Sources','Receiver','Background velocity')
title(sprintf('Spatial dof:%d, data dof:%d',n,source.n*receiver.n*freq.n))
axis([0 1 0 2])
drawnow

[m,out]=adjoint_state_1d(c,freq,source,receiver,maxit,noise,method);

figure
subplot(1,2,1)
% plot(x,c(x),'-b',x,sqrt(1./m(2:end-1)),'-r')
plot(x,1./c(x).^2,'-b',x,m(2:end-1),'-r')
legend('True squared slowness','Reconstruction')
axis([0 1 0 2])
subplot(1,2,2)
semilogy(out.J(out.J~=0))
title('Objective function')