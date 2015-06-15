function [m,out]=adjoint_state_2d(dom,all_freqs,sources,receivers,window_info,c_true,c_0,sigma,maxit)
%ADJOINT_STATE_2D   Adjoint-state method for full-waveform inversion in 2D.
%   [M,OUT] = ADJOINT_STATE_2D(DOM,ALL_FREQS,SOURCES,RECEIVERS,WINDOW_INFO,C_TRUE,C_0,SIGMA,MAXIT)
%   runs the adjoint-state implementation of full-waveform inversion. The
%   domain object DOM describes the physical domaind and grid used. The
%   vector ALL_FREQS is a list of all desired frequencies in hertz. The two
%   row matrices SOURCES and RECEIVERS contain x and y locations of sources
%   and receivers in the physical domain, and should be generated using
%   sources_and_receivers. Window_info is a struct containing information
%   about windowing; see dom.window for more. C_TRUE and C_0 are the true
%   and initial background velocities. SIGMA is the noise level and MAXIT
%   is the maximum number of permitted iterations.
%
%	Returned are a variable M, sometimes called the "squared slowness", and
%	a struct out with information coming directly from the LBFGS routine.

N = dom.N;									% total number of interior grid points
[win_inds,~,W] = dom.window(window_info);	% window indices and prolongation from window to domain
nf = length(all_freqs);						% number of frequencies
ns = size(sources,2);						% number of sources
nw = length(win_inds);						% number of pixels outside window
r_ind = dom.loc2ind(receivers);				% flat indices of receivers
I = speye(N);								% for padding/prolongation operators
E_rec = I(:,r_ind);							% padding operator
m_true = 1./c_true.^2;						% true sqaured slowness
m0 = m_true;								% initial squared slowness, true value inside window
m0(win_inds) = 1./c_0(win_inds).^2;			% set values inside of window to initial value
b = generate_sources(dom,sources);			% vectors corresponding to source locations
opt.Niter = maxit;							% setup opt for LBFGS

% Measure noisy data of true fields at receiver locations
d = generate_seismic_data(dom,sources,receivers,all_freqs,m_true,sigma);

% Run LBFGS with no initial perturbation
[dm,out] = lbfgs(@adjoint_state_gradient,zeros(nw,1),opt);

% Return updated m
m = m0+W*dm;

	function [J,DJ] = adjoint_state_gradient(dm)
		% Compute the gradient of J(m) = .5*|| S*u(m) - d ||^2 via the adjoint state method
		fprintf('Solving Helmholtz, frequency: ')
		DJ = zeros(nw,1);	% gradient
		J = 0;				% value of objective function
		% Parallelize over frequencies
		parfor ii = 1:nf
			fprintf('%d, ',ii)
			A = invertA(helmholtz_2d(m0+W*dm,all_freqs(ii),dom),1);
			for jj = 1:ns
				u = A.apply(b(:,jj)); %#ok<PFBNS>	% Forward solve
				bq = E_rec*(u(r_ind)-d(:,ii,jj));	% Padding
				q = A.applyt(bq);					% Adjoint solve
				% Update the gradient
				DJ = DJ+(2*pi*all_freqs(ii))^2*real(u(win_inds).*conj(q(win_inds)));
				% Update the objective value
				J = J+norm(u(r_ind)-d(:,ii,jj))^2;
			end
		end
		fprintf('\n')
	end
end