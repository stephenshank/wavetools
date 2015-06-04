% NOTICE: This script can be configured to save data from different runs.
% Simply set save_data = 1 and alter data_dir on line 12 below.

%% Clear and close everything, open parallel pool
clear all
close all
clc
if isempty(gcp('nocreate')), parpool; end

%% Key parameters (designed to be edited quickly)
window_type = 'all';									% change to get rectangular window around receivers
save_data = 0;											% determine whether we save info
data_dir = '~/Documents/MATLAB/wavetools/marm_data/';	% directory to save data in (if save_data=1)
fmin = 1;												% minimum frequency
fmax = 20;												% maximum frequency
nf = 12;												% number of frequencies
ctr = 1;												% contrast

%% Remaining parameters
nx = ceil(384*10*fmax/122);								% number of gridpoints in x direction
ny = ceil(10*fmax);										% number of gridpoints in y direction
wpml = .1;												% width of PML
freqs = linspace(fmin,fmax,nf);							% frequencies
ns = ceil(sqrt(nx*ny/(3*nf)));							% number of sources, in 1:3 ratio to receivers
nr = 3*ns;												% number of receivers
c_vec = [1 ctr+1];										% vector for colorbar
sigma = 1e-5;											% noise level
maxit = 5;												% maximum number of LBFGS iterations
smoothness = 10;										% intensity of smoothing filter
dom = domain([0 384/122 0 1],[nx ny]);					% rectangular domain
source_info.type = 'hline';								% create a horizontal line of sources
source_info.bounds = [.2 2.9];							% left and right endpoints of sources
source_info.height = .8;								% height of sources
sources = sources_and_receivers(ns,source_info);		% x and y locations of sources
receiver_info.type = 'hline';							% create a horizontal line of receivers
receiver_info.bounds = [.2 2.9];						% left and right endpoints of receivers
receiver_info.height = .75;								% height of receivers
receivers = sources_and_receivers(nr,receiver_info);	% x and y locations of receivers
if strcmp(window_type,'all')							% decide whether or not to window
	window_info.type = 'all';							% take all points
else
	window_info.type = 'rectangle_inner';				% create a rectangular window
	window_info.bounds = [.1 3 .7 .85];					% boundary of this window
end
[win_inds,W] = dom.window(window_info);					% indices of elements outside of window
pml_info.type = 'pml';									% pml info, for plotting purposes
pml_info.width = wpml;									% width of pml
[~,PML] = dom.window(pml_info);							% indicates whether a pixel is inside PML
c_true = dom.M2m(marmousi(dom,'hard',ctr));				% true background velocity
c0 = dom.M2m(smooth(dom.m2M(c_true),smoothness));		% initial estimate for background velocity

%% Saving data to disk
% Determine how many computational experiments have been run
if save_data
	exper_num = 1; %#ok<*UNRCH>
	while exist(sprintf('%s/setup%d.fig',data_dir,exper_num),'file')
		exper_num = exper_num+1;
	end
end

%% Preliminary plots of experimental setup
figure

% True background velocity, sources and receivers
subplot(2,4,1:2)
dom.imagesc(c_true);
hold on
dom.plot(receivers(1,:),receivers(2,:),'rx')
dom.plot(sources(1,:),sources(2,:),'go')
hold off
axis image
title('Problem setup')

% Initial background velocity
subplot(2,4,3:4)
dom.imagesc(c0)
axis image
title('Initial background velocity')
caxis(c_vec)

% Sources, receivers and PML
subplot(2,4,5:6)
imagesc(dom.m2M(c_true)+ctr*PML)
hold on
dom.plot(receivers(1,:),receivers(2,:),'rx')
dom.plot(sources(1,:),sources(2,:),'go')
hold off
title('With PML')
axis image
caxis([c_vec(1) c_vec(2)])

% Sources, receivers and window
subplot(2,4,7:8)
imagesc(dom.m2M(c_true)+ctr*flipud(W))
hold on
dom.plot(receivers(1,:),receivers(2,:),'rx')
dom.plot(sources(1,:),sources(2,:),'go')
hold off
title('With window')
axis image
caxis([c_vec(1) c_vec(2)])

% Format and save to disk
set(gcf,'Position',[11 340 1362 462])
if save_data
	savefig(sprintf('marm_data/setup%d.fig',exper_num))
	diary(sprintf('marm_data/output%d.txt',exper_num))
end

%% Reconstruction via adjoint state method
fprintf('Spatial dof:%d, inverse prob. dof:%d\n',dom.N,nf*ns*nr)
tic
[m,out] = adjoint_state_2d(dom,freqs,sources,receivers,window_info,c_true,c0,sigma,maxit);
cpu_time = toc;
fprintf('Running time... %.2f (s).\n',cpu_time)

%% Post computation reconstrution and objective function
figure 

% True background velocity
subplot(3,1,1)
dom.imagesc(c_true)
caxis(c_vec)
axis image
title('True model')
colorbar

% Initial approximation to background velocity
subplot(3,1,2)
dom.imagesc(c0)
caxis(c_vec)
axis image
title('Initial approximation')
colorbar

% Reconstruction via FWI
subplot(3,1,3)
if sum(m<0)>0, warning('Negative values encountered in m!'), end
dom.imagesc(sqrt(1./abs(m)))
caxis(c_vec)
axis image
title('Reconstruction')
colorbar
if save_data
	savefig(sprintf('marm_data/reconstruction%d.fig',exper_num))
end

figure
% Objective function
semilogy(out.J(out.J~=0))
axis square
title('Objective')
if save_data
	savefig(sprintf('marm_data/objective%d.fig',exper_num))
end

% Save variables for reproducibility
if save_data
	save(sprintf('marm_data/data%d',exper_num),'dom','ns','nr','freqs','source_info','receiver_info', ...
		'window_info','c_true','c0','ctr','sigma','maxit','m','out','cpu_time')
	diary off
end