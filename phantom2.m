function P=phantom2(n,p)

npad=ceil(p*n/2);
nint=n-2*npad;
P=zeros(n,n);
P(npad+1:n-npad,npad+1:n-npad)=phantom(nint);