function [K,M,omega] = helmholtz_2d(m,f,dom,pml)

xm = dom.xmin;				% left x boundary
xM = dom.xmax;				% right x boundary
ym = dom.ymin;				% bottom y boundary
yM = dom.ymax;				% top y boundary
nx = dom.nx;				% number of gridpoints in x direction
ny = dom.ny;				% number of gridpoints in y direction
hx = dom.hx;				% grid spacing in x direction
hy = dom.hy;				% grid spacing in y direction
x = dom.x;					% gridpoints in x direction
y = dom.y;					% gridpoints in y direction
omega = 2*pi*f;				% angular frequency
N = dom.N;					% total number of gridpoints
if nargin<4	|| isempty(pml)	% default PML parameters
	pml.width = .1;			% width
	pml.intensity = 1e4;	% intensity
end
xpml = pml.width;			% copy for readability
mpml = pml.intensity;		% copy for readability

% sigma functions for PML, as in 2008 Erlangga survey "Advances in..."
sigmax = @(x,y) mpml*((x<(xm+xpml)).*(x-(xm+xpml)).^2 + (x>(xM-xpml)).*(x-(xM-xpml)).^2);
sigmay = @(x,y) mpml*((y<(ym+xpml)).*(y-(ym+xpml)).^2 + (y>(yM-xpml)).*(y-(yM-xpml)).^2);

% corresponding alpha and beta functions
alpha = @(x,y) 1./(1+sigmax(x,y)/(1i*omega));
beta = @(x,y) 1./(1+sigmay(x,y)/(1i*omega));

xmids = ([xm;x]+[x;xM])/2;	% midpoints of x gridpoints, for alpha and beta terms
ymids = ([ym;y]+[y;yM])/2;	% midpoints of y gridpoints, for alpha and beta terms
Ax = cell(ny,1);			% diagonal blocks of A containing discretization of -d(alpha*du/dx)/dx
betas = zeros(N,3);			% diagonals of A containing discretization of -d(beta*du/dy)/dy

% Assemble contributions of modified Laplacian
for i = 1:ny
	% Discretization of -alpha*d(alpha*du/dx)/dx
	tridiag = [ -[alpha(xmids(2:end-1),y(i));0] alpha(xmids(1:end-1),y(i))+alpha(xmids(2:end),y(i)) ...
		-[0;alpha(xmids(2:end-1),y(i))] ];
	% Next two lines are nontrivial... need to consider which elements spdiags discards
	alpha_mult = [ [alpha(x(2:end),y(i));0] alpha(x,y(i)) [0;alpha(x(1:end-1),y(i))] ];
	Ax{i} = spdiags(alpha_mult.*tridiag,-1:1,nx,nx);
	
	% Discretization of -beta*d(beta*du/dy)/dy
	ind = 1+(i-1)*nx:i*nx;
	if i<ny
		betas(ind,1) = -beta(x(1),y(i+1))*beta(x,ymids(i+1));
	end
	betas(ind,2) = beta(x(1),y(i))*(beta(x,ymids(i))+beta(x,ymids(i+1)));
	if i>1
		betas(ind,3) = -beta(x(1),y(i-1))*beta(x,ymids(i));
	end
end

if nargout==1
	K = blkdiag(Ax{:})/hx^2 + ...					% x derivatives in modified Laplacian
		spdiags(betas,[-nx 0 nx],N,N)/hy^2 - ...	% y derivatives in modified Laplacian
		omega^2*spdiags(m,0,N,N);					% discretization of m*omega^2*s1*s2
else
	K = blkdiag(Ax{:})/hx^2 + ...			% x derivatives in modified Laplacian
		spdiags(betas,[-nx 0 nx],N,N)/hy^2;	% y derivatives in modified Laplacian
	M =	spdiags(m,0,N,N);			% discretization of m*omega^2*s1*s2
end