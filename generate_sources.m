function b = generate_sources(dom,sources)

ns = size(sources,2);		% number of sources
b = zeros(dom.N,ns);		% initialize right-hand sides
% Create discrete deltas at source location
for j = 1:ns
	b(:,j) = dom.pt_src(sources(1,j),sources(2,j));
end