function M=marmousi(dom,type,ctr,smoothness)
% smoothness is an optional parameter, determines smoothing amount of Gaussian filter

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
[Xq,Yq]=ndgrid(linspace(1,122,dom.ny),linspace(1,384,dom.nx));
Mgi=griddedInterpolant(Xg,Yg,M,'nearest');
M=Mgi(Xq,Yq);
if nargin==4 && smoothness>0
	M=smooth(M,smoothness);
end
M=1+ctr*(M-min(M(:)))/(max(M(:))-min(M(:)));