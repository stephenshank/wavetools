function [m,out]=adjoint_state_2d(dom,freqs,sources,receivers,window_info,ctrue,c0,sigma,maxit)
%ADJOINT_STATE_2D   Adjoint-state method for full-waveform inversion in 2D.

if isempty(gcp('nocreate')), parpool; end

N = dom.N;								% total number of interior grid points
[win_inds,~,W] = dom.window(window_info);
nf = length(freqs);						% number of frequencies
ns = size(sources,2);					% number of sources
nr = size(receivers,2);					% number of receivers
nw = length(win_inds);					% number of pixels outside window
r_ind = dom.loc2ind(receivers);			% flat indices of receivers
b = zeros(N,ns);						% initialize right-hand sides
d = zeros(nr,nf,ns);					% seismic data
I = speye(N);							% for padding/prolongation operators
E_rec = I(:,r_ind);						% padding operator
mtrue = 1./ctrue.^2;					% true sqaured slowness
m0 = mtrue;								% initial squared slowness, true value inside window
m0(win_inds) = 1./c0(win_inds).^2;		% set values outside of window to initial value

% Create discrete deltas at source location
for j = 1:ns
	b(:,j) = dom.pt_src(sources(1,j),sources(2,j));
end

% Generate true data
fprintf('Generating data... frequency: ')
rng(1)
parfor i = 1:nf
	fprintf('%d, ',i)
	A_true = invertA(helmholtz_2d(mtrue,freqs(i),dom));
	for j = 1:ns
		u_true = A_true.apply(b(:,j)); %#ok<*PFBNS>
		d(:,i,j) = u_true(r_ind)+sigma*(rand(nr,1)+1i*rand(nr,1));
	end
end

fprintf('\n')
clear u_true A_true
opt.Niter = maxit;
[dm,out] = lbfgs(@adjoint_state_gradient,zeros(nw,1),opt);
m = m0+W*dm;

	function [J,DJ] = adjoint_state_gradient(dm)
		% Compute the gradient of J(m) = .5*|| S*u(m) - d ||^2 via the adjoint state method
		fprintf('Solving Helmholtz, frequency: ')
		DJ = zeros(nw,1);	% gradient
		J = 0;				% value of objective function
		parfor ii = 1:nf
			fprintf('%d, ',ii)
			A = invertA(helmholtz_2d(m0+W*dm,freqs(ii),dom),1);
			for jj = 1:ns
				u = A.apply(b(:,jj));
				bq = E_rec*(u(r_ind)-d(:,ii,jj));
				q = A.applyt(bq);
				DJ = DJ+(2*pi*freqs(ii))^2*real(u(win_inds).*conj(q(win_inds)));
				J = J+norm(u(r_ind)-d(:,ii,jj))^2;
			end
		end
		fprintf('\n')
	end
end