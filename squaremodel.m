function P=squaremodel(n,ctr,bounds,smoothness)

if nargin<2, ctr=1; end
if nargin<3, bounds=[.3 .7 .3 .7]; end
[X,Y] = meshgrid(linspace(0,1,n));
P = 1 + ctr*( bounds(1) < X & X < bounds(2) & bounds(3) < Y & Y < bounds(4));
if nargin==4 && smoothness>0
	P = smooth(P,smoothness);
end