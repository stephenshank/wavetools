clear all
close all 
clc
% function adjoint_state_1d
fmin=1e-4;						% minimum frequency
fmax=20;						% maximum frequency
maxit=100;
err=zeros(maxit,1);
n=30*fmax;						% number of interior grid points
h=1/(n+1);						% spatial resolution
nx=n+2;							% total number of grid points
nf=20;							% number of frequencies
ns=10;							% number of sources
nr=8;							% number of receivers
sigma=1e-8;						% noise level
freqs=linspace(fmin,fmax,nf);	% vector of all frequencies
x=linspace(0,1,nx);				% grid points

% Locations of sources and receivers
sl=[.05 .3];
rl=[.1 .35];
sr=[.7 .95];
rr=[.65 .9];
xs_ind=[round(linspace(sl(1),sl(2),ns/2)/h) round(linspace(sr(1),sr(2),ns/2)/h)];
xr_ind=[round(linspace(rl(1),rl(2),nr/2)/h) round(linspace(rr(1),rr(2),nr/2)/h)];

% sint=[.1 .3];
% rint=[.65 .9];
% xs_ind=round(linspace(sint(1),sint(2),ns)/h);
% xr_ind=round(linspace(rint(1),rint(2),nr)/h);

c=@(x) 1+.2*exp(-100*(x-.5).^2);	% true background velocity
% c=@(x) 1 + .1*((x>.4) & (x<.6));
% c=@(x) 1 + .2*(((x>.4) & (x<.45)) | ((x>.5) & (x<.55)));
m=ones(nx,1);					% initial guess for background squared slowness
mtrue=(1./c(x).^2)';			% true background squared slowness
% mtrue=m;
% mtrue=m;			% true background squared slowness
u=zeros(nx,1);				% field
q=zeros(nx,1);				% adjoint state
b=zeros(nx,ns);					% right-hand sides
d=zeros(nr,nf,ns);				% data
I=speye(nx);
E=I(:,xr_ind);

% Sanity checks
% assert(ns==length(xs_ind))
% assert(nr==length(xr_ind))

% Plot background velocity, along with locations of sources and receivers
figure
plot(x(xs_ind),ones(size(xs_ind)),'rx',x(xr_ind),ones(size(xr_ind)),'bo',x,sqrt(1./mtrue),'-k')
legend('Sources','Receiver','Background velocity')
title(sprintf('Spatial dof:%d, data dof:%d',n,ns*nr*nf))
axis([0 1 0 2])
pause

% create right-hand sides and data
rng(1)
for i=1:nf
	A=helmholtz_fos_1d(mtrue,freqs(i));
	for j=1:ns
		b(xs_ind(j),j)=1/h;
		u=A\b(:,j);
		d(:,i,j)=u(xr_ind)+sigma*(rand+1i*rand);
	end
end

figure
% Gradient descent loop
for k=1:maxit
% 	plot(x(2:end-1),sqrt(1./m(2:end-1)),'-rx',x,sqrt(1./mtrue),'-b')
	if k>1
		plot(x,m,'-rx',x,mtrue,'-b')
		title(sprintf('Iter:%d, misfit:%.3e',k,err(k-1)))
		axis([0 1 -2 2])
		drawnow
	end
% 	pause
	% Compute gradient via adjoint state method
	f=zeros(nx,1);
	for i=1:nf
		A=helmholtz_fos_1d(m,freqs(i));
		for j=1:ns
			u=A\b(:,j);
			bq=E*(u(xr_ind)-d(:,i,j));
% 			q=A\(E*conj(u(xr_ind)-d(:,i,j)));
			q=(A')\bq;
% 			fprintf('|<A^*q,b>-<q,Ab>|=%.5e\n',abs((A'*q)'*bq-q'*(A*bq)))
% 			plot(x,real(u),'-b',x,real(q),'-r',x(xs_ind),zeros(size(xs_ind)),'gx',x(xr_ind),zeros(size(xr_ind)),'ko')
% 			title(sprintf('freq:%d, source:%d',i,j))
% 			pause
			f=f+(2*pi*freqs(i))^2*real(u.*conj(q));
			err(k)=err(k)+norm(u(xr_ind)-d(:,i,j))^2;
		end
		if any(isnan(m)), error('NaNs in m!'), end
	end
	disp(norm(f))
	alpha=(f'*m)/(f'*f);
	m=m-alpha*f;
end