function P=phantom2(n,p,ctr,smoothness)

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