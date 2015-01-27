function A = helmholtz_2d(m,f,n,pml)
omega = 2*pi*f;				% angular frequency
N = n^2;					% total number of gridpoints
if nargin<4					% default PML parameters
	pml.width = .1;			% width
	pml.intensity = 1e4;	% intensity
end
h = 1/(n+1);				% grid spacing
x = (h:h:1-h)';				% vector of grid points
xpml = pml.width;			% copy for readability
mpml = pml.intensity;		% copy for readability

% sigma functions for PML, as in 2008 Erlangga survey "Advances in..."
sigmax = @(x,y) mpml*((x<xpml).*(x-xpml).^2 + (x>(1-xpml)).*(x-(1-xpml)).^2);
sigmay = @(x,y) mpml*((y<xpml).*(y-xpml).^2 + (y>(1-xpml)).*(y-(1-xpml)).^2);

% corresponding s functions
sx = @(x,y) 1+sigmax(x,y)/(1i*omega);
sy = @(x,y) 1+sigmay(x,y)/(1i*omega);

% Due to PML, Laplacian becomes - d(alpha*du/dx)/dx - d(beta*du/dy)/dy 
alpha = @(x,y) sy(x,y)./sx(x,y);
beta = @(x,y) sx(x,y)./sy(x,y);

[X,Y] = meshgrid(x);	% grid points in meshgrid format
xm = ([0;x]+[x;1])/2;	% midpoints of gridpoints, for alpha and beta terms
Ax = cell(n,1);			% diagonal blocks of A containing discretizaiton of -d(alpha*du/dx)/dx
betas = zeros(N,3);		% diagonals of A containing discretizaiton of -d(beta*du/dy)/dy

% Assemble contributions of modified Laplacian
for i = 1:n
	ind = 1+(i-1)*n:i*n;
	if i<n
		betas(ind,1) = -beta(x,xm(i+1));
	end
	betas(ind,2) = beta(x,xm(i))+beta(x,xm(i+1));
	if i>1
		betas(ind,3) = -beta(x,xm(i));
	end
	tridiag = [ -[alpha(xm(2:end-1),x(i));0] alpha(xm(1:end-1),x(i))+alpha(xm(2:end),x(i)) ...
		-[0;alpha(xm(2:end-1),x(i))] ];
	Ax{i} = spdiags(tridiag,-1:1,n,n);
end
A = (spdiags(betas,[-n 0 n],N,N) + blkdiag(Ax{:}))/h^2 - ...		% modified Laplacian
	omega^2*spdiags(vec((sx(X,Y).*sy(X,Y)).'.*flipud(m)'),0,N,N);	% discretization of m*omega^2*s1*s2