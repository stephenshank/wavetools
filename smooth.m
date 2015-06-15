function S = smooth(M,smoothness)
% SMOOTH   Apply a smoothing filter.
%   S = SMOOTH(M,SMOOTHNESS) applies a Gaussian smoothing filter to the
%   model represented by the matrix M with smoothing intensity given by
%   SMOOTHNESS. The filter is 50 pixels by 50 pixels.

S=imfilter(M,fspecial('gaussian',[50 50],smoothness),'replicate');