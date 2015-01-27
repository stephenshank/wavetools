clear all
close all
clc

tic

freq.min=1;				% minimum frequency
freq.max=5;				% maximum frequency
freq.n=10;				% number of frequencies
ns=25;					% number of sources
nr=35;					% number of receivers
n=ceil(10*freq.max);	% number of grid points (10 per wavelength)
ctr=.3;					% contrast
sigma=1e-6;				% noise level
maxit=50;				% maximum number of iterations

% Shepp-Logan
c_true=@(n) phantom2(n,.3,ctr);
c0=@(n) ones(n);
source_data.type='circle';
source_data.center=[.5 .5];
source_data.radius=.35;
sources=sources_2d(ns,source_data);
receiver_data.type='circle';
receiver_data.center=[.5 .5];
receiver_data.radius=.4;
receivers=sources_2d(nr,receiver_data);

% Marmousi
% c_true=@(n) marmousi(n,'hard',ctr);
% c0=@(n) marmousi(n,'smooth',ctr);
% source_data.type='hline';
% source_data.left=.15;
% source_data.right=.85;
% source_data.height=.85;
% sources=sources_2d(ns,source_data);
% receiver_data1.type='vline';
% receiver_data1.top=.1;
% receiver_data1.bottom=.9;
% receiver_data1.pos=.1;
% receiver_data2.type='vline';
% receiver_data2.top=.1;
% receiver_data2.bottom=.9;
% receiver_data2.pos=.9;
% receivers=[sources_2d(nr,receiver_data1);sources_2d(nr,receiver_data2)];

fprintf('Spatial dof:%d, inverse prob. dof:%d\n',n^2,freq.n*ns*nr)
[m,out]=adjoint_state_2d(freq,sources,receivers,c_true,maxit,c0,sigma);

f=figure;

subplot(2,2,1)
imagesc(c_true(n));
c_vec=caxis;
hold on
s_xind=loc2ind(n,sources(:,1)); s_yind=loc2ind(n,sources(:,2));
r_xind=loc2ind(n,receivers(:,1)); r_yind=loc2ind(n,receivers(:,2));
plot(r_xind,n-r_yind+1,'rx',s_xind,n-s_yind+1,'go')
axis square
legend('Receivers','Sources')
title('Problem setup')

subplot(2,2,2)
imagesc(c0(n))
axis square
title('Initial background velocity')
caxis(c_vec)

subplot(2,2,3)
if sum(m<0)>0, warning('Negative values encountered in m!'), end
imagesc(sqrt(1./flipud(reshape(abs(m),n,n)')))
caxis(c_vec)
axis square
title('Reconstruction')

subplot(2,2,4)
semilogy(out.J(out.J~=0))
axis square
title('Objective')

toc