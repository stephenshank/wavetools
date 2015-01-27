function P=phantom2(n,p,ctr)

if nargin<2, p=0; end
if nargin<3, ctr=1; end
npad=ceil(p*n/2);
nint=n-2*npad;
P=ones(n,n);
P(npad+1:n-npad,npad+1:n-npad)=P(npad+1:n-npad,npad+1:n-npad)+ctr*phantom(nint);