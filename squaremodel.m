function P = squaremodel(dom,ctr,bounds,smoothness)
% SQUAREMODEL Generate a model which is a piecewise constant square.
%   P = SQUAREMODEL(DOM,CTR,BOUNDS,SMOOTHNESS) returns a matrix P which
%   represents a velocity model of a square on the grid associated with the
%   domain object DOMAIN. The remaining parameters are optional.
%
%   CTR is the constrast of the resulting model. SMOOTHNESS will
%   will apply a Gaussian smoothing filter of 50 pixels for a given
%   smoothing insensity. BOUNDS is a vector of the form
%   [xmin xmax ymin ymax] which give the bounds of the square, which is set
%   to [.3 .7 .3 .7] by default.
%
%   See also DOMAIN, SMOOTH.

n = dom.nx;
if nargin<2, ctr=1; end
if nargin<3, bounds=[.3 .7 .3 .7]; end
[X,Y] = meshgrid(linspace(0,1,n));
P = 1 + ctr*( bounds(1) < X & X < bounds(2) & bounds(3) < Y & Y < bounds(4));
if nargin==4 && smoothness>0
	P = smooth(P,smoothness);
end