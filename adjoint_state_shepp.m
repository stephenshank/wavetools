%% Clear and close everything, open parallel pool
clear all
close all
clc
if isempty(gcp('nocreate')), parpool; end

%% Problem setup and parameters
fmin=1;							% minimum frequency
fmax=10;						% maximum frequency
nf=10;							% number of frequencies
wpml=.1;						% width of PML
freqs=linspace(fmin,fmax,nf);	% frequencies
ns=20;							% number of sources
nr=50;							% number of receivers
ctr=.5;							% contrast
sigma=1e-5;						% noise level
maxit=25;						% maximum number of LBFGS iterations
smoothness=.5;					% intensity of smoothing filter

%% Shepp-Logan phantom and source/receiver/window configuration
nx=ceil(10*fmax);				% number of gridpoints in a given direction
dom=domain([0 1 0 1],[nx nx]);	% square domain with equispaced gridpoints
source_info.type='ring';		% create a ring of sources
source_info.center=[.5 .5];		% center of ring of sources
source_info.radius=.35;			% radius of ring of sources
sources=sources_and_receivers(ns,source_info);
receiver_info.type='ring';		% create a ring of receivers
receiver_info.center=[.5 .5];	% center of ring of receivers
receiver_info.radius=.4;		% radius of ring of receivers
receivers=sources_and_receivers(nr,receiver_info);
window_info.type='disk';
window_info.center=[.5 .5];
window_info.radius=.3;

[~,W]=dom.window(window_info);
c_true=dom.M2m(phantom2(nx,.35,ctr,smoothness));
c0=dom.M2m(ones(nx));
pml_info.type='pml';
pml_info.width=wpml;
[~,PML]=dom.window(pml_info);

%% Preliminary plots of experimental setup
figure

subplot(2,2,1)
dom.imagesc(c_true);
c_vec=caxis;
hold on
dom.plot(receivers(1,:),receivers(2,:),'rx')
dom.plot(sources(1,:),sources(2,:),'go')
hold off
axis square
title('Problem setup')

subplot(2,2,2)
dom.imagesc(c0)
axis square
title('Initial background velocity')
caxis(c_vec)

subplot(2,2,3)
imagesc(dom.m2M(c_true)+PML)
hold on
dom.plot(receivers(1,:),receivers(2,:),'rx')
dom.plot(sources(1,:),sources(2,:),'go')
hold off
title('With PML')
axis square
caxis([c_vec(1) c_vec(2)])

subplot(2,2,4)
imagesc(dom.m2M(c_true)+flipud(W))
hold on
dom.plot(receivers(1,:),receivers(2,:),'rx')
dom.plot(sources(1,:),sources(2,:),'go')
hold off
title('With window')
axis square
caxis([c_vec(1) c_vec(2)+1])

drawnow

%% Reconstruction via adjoint state method
fprintf('Spatial dof:%d, inverse prob. dof:%d\n',dom.N,nf*ns*nr)
tic
[m,out]=adjoint_state_2d(dom,freqs,sources,receivers,window_info,c_true,c0,sigma,maxit);
toc

%% Post computation reconstrution and objective function
figure 

subplot(1,2,1)
if sum(m<0)>0, warning('Negative values encountered in m!'), end
dom.imagesc(sqrt(1./m))
caxis(c_vec)
axis square
title('Reconstruction')

subplot(1,2,2)
semilogy(out.J(out.J~=0))
axis square
title('Objective')