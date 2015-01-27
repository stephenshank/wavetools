function adjoint_state_2d
clear all
close all
clc
fmin=1;							% minimum frequency
fmax=15;						% maximum frequency
n=10*fmax;						% number of interior grid points in one direction
N=n^2;							% total number of interior grid points
h=1/(n+1);						% spatial resolution
nf=10;							% number of frequencies
ns=50;							% number of sources
nr=50;							% number of receivers
sigma=1e-2;						% noise level
freqs=linspace(fmin,fmax,nf);	% vector of all frequencies
% x=h:h:1-h;						% grid points
ctr=1;							% contrast
c_true = 1+ctr*phantom2(n,.3);

% Locations of sources and receivers
s_center=[.5 .5];
s_radius=.35;
s_theta=2*pi*(0:ns-1)/ns;
s_xloc=s_center(1)+s_radius*cos(s_theta);
s_yloc=s_center(2)+s_radius*sin(s_theta);
s_xind=round(s_xloc/h);
s_yind=round(s_yloc/h);
s_fi=sub2ind([n n],s_xind,s_yind);

r_center=[.5 .5];
r_radius=.4;
r_theta=2*pi*(0:ns-1)/ns;
r_xloc=r_center(1)+r_radius*cos(r_theta);
r_yloc=r_center(2)+r_radius*sin(r_theta);
r_xind=round(r_xloc/h);
r_yind=round(r_yloc/h);
r_fi=sub2ind([n n],r_xind,r_yind);

figure
imagesc(c_true)
c_vec=caxis;
hold on
plot(r_xind,n-r_yind+1,'rx',s_xind,n-s_yind+1,'go')
axis square
legend('Receivers','Sources')
title(sprintf('Spatial dof:%d, inverse prob. dof:%d',N,nf*ns*nr))
drawnow

m=ones(N,1);					% initial guess for background squared slowness
b=zeros(N,ns);					% right-hand sides
d=zeros(nr,nf,ns);				% data
I=speye(N);
E=I(:,r_fi);

% Create right-hand sides and data
fprintf('Generating data... frequency: ')
rng(1)
for i=1:nf
	fprintf('%d, ',i)
	A_true=invertA(helmholtz_2d(1./c_true.^2,freqs(i),n));
	for j=1:ns
		if i==1, b(s_fi(j),j)=1/h^2; end
		u_true=A_true.apply(b(:,j));
		d(:,i,j)=u_true(r_fi)+sigma*(rand+1i*rand);
	end
end
fprintf('\n')
clear u_true A_true
opt.Niter=40;
[m,out]=lbfgs(@adjoint_state_gradient,m,opt);

figure
subplot(1,2,1)
imagesc(sqrt(1./flipud(reshape(m,n,n)')))
caxis(c_vec)
title('Reconstruction')
subplot(1,2,2)
semilogy(out.J(out.J~=0))
title('Objective')

	function [J,DJ]=adjoint_state_gradient(m)
		fprintf('Solving Helmholtz, frequency: ')
		J=0;
		DJ=zeros(N,1);
		for ii=1:nf
			fprintf('%d, ',ii)
			A=invertA(helmholtz_2d(flipud(reshape(m,n,n)'),freqs(ii),n),1);
			for jj=1:ns
				u=A.apply(b(:,jj));
				bq=E*(u(r_fi)-d(:,ii,jj));
				q=A.applyt(bq);
				DJ=DJ+(2*pi*freqs(ii))^2*real(u.*conj(q));
				J=J+norm(u(r_fi)-d(:,ii,jj))^2;
			end
		end
		fprintf('\n')
	end
end