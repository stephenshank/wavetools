function b = pt_src(n,x,y)

h = 1/(n+1);				% mesh width
xs = round(x/h);			% 1d index of x coord. of point source
ys = round(y/h);			% 1d index of y coord. of point source
is = sub2ind([n n],xs,ys);	% global index of point source in reshaped vector
b = zeros(n^2,1);			% initialize source vector
b(is) = 1/h^2;				% discrete Dirac delta