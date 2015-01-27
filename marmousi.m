function M=marmousi(n,type,ctr)

if nargin<2, type='hard'; end
if nargin<3, ctr=1; end
switch type
	case 'hard'
		load marmhard.dat
		M=reshape(marmhard,122,[]);
	case 'smooth'
		load marmsmooth.dat
		M=reshape(marmsmooth,122,[]);
end
[Xg,Yg]=ndgrid(1:122,1:384);
[Xq,Yq]=ndgrid(linspace(1,122,n),linspace(1,384,n));
Mgi=griddedInterpolant(Xg,Yg,M);
M=Mgi(Xq,Yq);
M=M/min(M(:));
M=1+ctr*(M-1);