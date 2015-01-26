function M=marmousi(n)

load marmhard.dat
M=reshape(marmhard,122,[]);
[Xg,Yg]=ndgrid(1:122,1:384);
[Xq,Yq]=ndgrid(linspace(1,122,n),linspace(1,384,n));
Mgi=griddedInterpolant(Xg,Yg,M);
M=Mgi(Xq,Yq);
M=M/min(M(:));