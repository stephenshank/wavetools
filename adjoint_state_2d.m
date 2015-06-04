function [m,out]=adjoint_state_2d(dom,all_freqs,sources,receivers,window_info,ctrue,c0,sigma,maxit)
%ADJOINT_STATE_2D   Adjoint-state method for full-waveform inversion in 2D.

N = dom.N;									% total number of interior grid points
[win_inds,~,W] = dom.window(window_info);
nf = length(all_freqs);						% number of frequencies
ns = size(sources,2);						% number of sources
nw = length(win_inds);						% number of pixels outside window
r_ind = dom.loc2ind(receivers);				% flat indices of receivers
I = speye(N);								% for padding/prolongation operators
E_rec = I(:,r_ind);							% padding operator
m_true = 1./ctrue.^2;						% true sqaured slowness
m0 = m_true;								% initial squared slowness, true value inside window
m0(win_inds) = 1./c0(win_inds).^2;			% set values outside of window to initial value
b = generate_sources(dom,sources);
d = generate_seismic_data(dom,sources,receivers,all_freqs,m_true,sigma);

fprintf('\n')
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
			A = invertA(helmholtz_2d(m0+W*dm,all_freqs(ii),dom),1);
			for jj = 1:ns
				u = A.apply(b(:,jj)); %#ok<PFBNS>
				bq = E_rec*(u(r_ind)-d(:,ii,jj));
				q = A.applyt(bq);
				DJ = DJ+(2*pi*all_freqs(ii))^2*real(u(win_inds).*conj(q(win_inds)));
				J = J+norm(u(r_ind)-d(:,ii,jj))^2;
			end
		end
		fprintf('\n')
	end
end