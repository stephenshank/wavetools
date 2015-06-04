function d = generate_seismic_data(dom,sources,receivers,all_freqs,m_true,sigma)
nf = length(all_freqs);			% number of frequencies
ns = size(sources,2);			% number of sources
nr = size(receivers,2);			% number of receivers
d = zeros(nr,nf,ns);			% seismic data
r_ind = dom.loc2ind(receivers);	% flat indices of receivers

% Generate true data
fprintf('Generating data... frequency: ')
rng(1)	% set random seed
parfor i = 1:nf
	fprintf('%d, ',i)
	% Invert Helmholtz operator corresponding to true model
	A_true = invertA(helmholtz_2d(m_true,all_freqs(i),dom));
	for j = 1:ns
		% Apply to jth source
		u_true = A_true.apply(dom.pt_src(sources(1,j),sources(2,j))); %#ok<*PFBNS>
		% Gather seismic data
		d(:,i,j) = u_true(r_ind)+sigma*(rand(nr,1)+1i*rand(nr,1));
	end
end