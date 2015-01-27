function [m,misfit]=adjoint_state_2d

fmin=1;							% minimum frequency
fmax=20;						% maximum frequency
n=10*fmax;						% number of interior grid points in one direction
N=n^2;							% total number of interior grid points
h=1/(n+1);						% spatial resolution
nf=20;							% number of frequencies
ns=50;							% number of sources
nr=50;							% number of receivers
sigma=1e-2;						% noise level
freqs=linspace(fmin,fmax,nf);	% vector of all frequencies
x=h:h:1-h;						% grid points
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
imagesc()
hold on
plot(r_xind,n-r_yind+1,'rx',s_xind,n-s_yind+1,'go')
axis equal
legend('Receivers','Sources')

m=ones(N,1);					% initial guess for background squared slowness

b=zeros(N,ns);					% right-hand sides
d=zeros(nr,nf,ns);				% data
I=speye(N);
E=I(:,r_fi);

% Create right-hand sides and data
rng(1)
for i=1:nf
	A_true=helmholtz_fos_1d(1./c_true.^2,freqs(i));
	for j=1:ns
		b(s_fi(j),j)=1/h^2;
		u_true=A_true\b(:,j);
		d(:,i,j)=u_true(r_fi)+sigma*(rand+1i*rand);
	end
end
clear u_true A_true
% 
% [m,out]=LBFGS(@adjoint_state_gradient,m);
% misfit=out.J(out.J~=0);
% 
% 	function [J,DJ]=adjoint_state_gradient(m)
% 		J=0;
% 		DJ=zeros(nx,1);
% 		for ii=1:nf
% 			A=helmholtz_fos_1d(m,freqs(ii));
% 			for jj=1:ns
% 				u=A\b(:,jj);
% 				bq=E*(u(xr_ind)-d(:,ii,jj));
% 				q=(A')\bq;
% 				DJ=DJ+(2*pi*freqs(ii))^2*real(u.*conj(q));
% 				J=J+norm(u(xr_ind)-d(:,ii,jj))^2;
% 			end
% 		end
% 	end
% end