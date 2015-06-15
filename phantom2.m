function P = phantom2(dom,p,ctr,smoothness)
% PHANTHOM2   Generate Shepp-Logan Phantom with modifications.
%   P = PHANTOM2(N,P,CTR,SMOOTHNESS) returns a matrix P which represents
%   the Shepp-Logan model on the grid associated with the domain DOM. The
%   remaining parameters are optional.
%
%   P represents a percentage shrinkage of the Shepp-Logan phantom; the
%   resulting phantom will be shrunk down in order to fit within windows
%   and PMLs. CTR is the constrast of the resulting model. SMOOTHNESS will
%   will apply a Gaussian smoothing filter of 50 pixels for a given
%   smoothing insensity.
%
%   See also DOMAIN, SMOOTH.

n = dom.nx;
if nargin<2, p=0; end
if nargin<3, ctr=1; end
npad=ceil(p*n/2);
nint=n-2*npad;
P=zeros(n,n);
P(npad+1:n-npad,npad+1:n-npad)=phantom(nint);
if nargin==4 && smoothness > 0
	P=smooth(P,smoothness);
end
P=1+ctr*(P-min(P(:)))/(max(P(:))-min(P(:)));