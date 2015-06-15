function locs = sources_and_receivers(num,info)
% SOURCES_AND_RECEIVERS   Create locations of sources and receivers.
%   LOCS = SOURCES_AND_RECEIVERS(NUM,INFO) returns locations of
%   sources/receivers. The amount returned is dictacted by NUM, and their
%   arrangement is determined by the struct INFO.
%
%   INFO should also have a field type, which is a string that can either 
%   be ring, hline, or vline. If type is set equal to ring, INFO should
%   have fields center and radius that describe the associated ring. If
%   type is set to hline, INFO should must a field bounds which a is vector
%   describing left and right bounds of a horizontal line and a field
%   height which gives the of this line. Similary, if type is vline, INFO
%   must have a field BOUNDS which is a vector that gives the lower and
%   upper points on the associated vertical line, as well as a field POS
%   which gives the horizontal position of this line.

locs=zeros(2,num);
switch info.type
	case 'ring'
		theta=2*pi*(0:num-1)/num;
		locs(1,:)=info.center(1)+info.radius*cos(theta);
		locs(2,:)=info.center(2)+info.radius*sin(theta);
	case 'hline'
		locs(1,:)=linspace(info.bounds(1),info.bounds(2),num);
		locs(2,:)=info.height;
	case 'vline'
		locs(1,:)=info.pos;
		locs(2,:)=linspace(info.bounds(1),info.bounds(2),num);
end