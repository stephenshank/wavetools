function [m,out]=adjoint_state_2d(freq,sources,receivers,ctrue,maxit,c0,sigma)

if isempty(gcp('nocreate')), parpool; end
fmin=freq.min;							% minimum frequency
fmax=freq.max;						% maximum frequency
n=10*fmax;						% number of interior grid points in one direction
N=n^2;							% total number of interior grid points
h=1/(n+1);						% spatial resolution
nf=freq.n;							% number of frequencies
ns=size(sources,1);							% number of sources
nr=size(receivers,1);							% number of receivers
freqs=linspace(fmin,fmax,nf);	% vector of all frequencies

% Flat indices of sources and receivers
s_xind=loc2ind(n,sources(:,1));
s_yind=loc2ind(n,sources(:,2));
s_fi=sub2ind([n n],s_xind,s_yind);

r_xind=loc2ind(n,receivers(:,1));
r_yind=loc2ind(n,receivers(:,2));
r_fi=sub2ind([n n],r_xind,r_yind);

b=zeros(N,ns);					% right-hand sides
d=zeros(nr,nf,ns);				% data
I=speye(N);
E=I(:,r_fi);

% Create right-hand sides and data
for j=1:ns
	b(s_fi(j),j)=1/h^2;
end
fprintf('Generating data... frequency: ')
rng(1)
parfor i=1:nf
	fprintf('%d, ',i)
	A_true=invertA(helmholtz_2d(1./ctrue(n).^2,freqs(i),n));
	for j=1:ns
		u_true=A_true.apply(b(:,j)); %#ok<*PFBNS>
		d(:,i,j)=u_true(r_fi)+sigma*(rand+1i*rand);
	end
end
fprintf('\n')
clear u_true A_true
opt.Niter=maxit;
[m,out]=lbfgs(@adjoint_state_gradient,vec(flipud((1./c0(n).^2)')),opt);

	function [J,DJ]=adjoint_state_gradient(m)
		fprintf('Solving Helmholtz, frequency: ')
		J=0;
		DJ=zeros(N,1);
		parfor ii=1:nf
			fprintf('%d, ',ii)

			A=invertA(helmholtz_2d(flipud(reshape(m,n,n)'),freqs(ii),n));
			for jj=1:ns
				u=A.apply(b(:,jj));
				bq=E*(u(r_fi)-d(:,ii,jj));
				qbar=A.apply(conj(bq));
				DJ=DJ+(2*pi*freqs(ii))^2*real(u.*qbar);
				J=J+norm(u(r_fi)-d(:,ii,jj))^2;
			end
		end
		fprintf('\n')
	end
end