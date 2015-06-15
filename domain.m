classdef domain
%DOMAIN   Physical domain, discrete grid and related quantities.
%   DOM = DOMAIN(SPATIAL_DATA,SIZE_DATA) returns an object which has useful
%   methods for working PDEs on a two-dimensional domain that are
%   discretized on a grid in the interior of that domain.
%
%   SPATIAL_DATA is a vector with four entries [xmin,xmax,ymin,ymax] which
%   give the dimensions of the rectangular domain. SIZE_DATA is a vector
%   with two entries which gives the number of interior gridpoints in the
%   x- and y-directions, respectively.
%
%   domain Properties:
%      xmin - Minimal x-value describing rectangular domain.
%      xmax - Maximal x-value describing rectangular domain.
%      ymin - Minimal y-value describing rectangular domain.
%      ymax - Maximal y-value describing rectangular domain.
%      nx - Number of interior gridpoints in x-direction.
%      ny - Number of interior gridpoints in y-direction.
%      N - Total number of interior gridpoints.
%      hx - Grid spacing in x-direction.
%      hy - Grid spacing in y-direction.
%      x - Vector of x-values of interior gridpoints.
%      y - Vector of y-values of interior gridpoints.
%      X - Matrix of x-values of interior gridpoints (from meshgrid).
%      Y - Matrix of y-values of interior gridpoints (from meshgrid).
%
%   domain Methods:
%      vec2mat - Convert a vector to matrix on the grid.
%      mat2vec - Convert a matrix on the grid to a vector.
%      loc2sub - Convert locations of points on the grid to subscripts on
%         grid.
%      sub2ind - Convert subscripts on the grid to flat indices.
%      loc2ind - Convert locations of points on the grid to flat indices.
%      pt_src - Create a point source (Dirac delta) at a given location.
%      window - Gives relevant quantities for working with a window on the
%         given domain.
%      imagesc - Analogue of MATLAB's imagesc for working on the grid.
%      plot - Analogue of MATLAB's plot for working on the grid.
%      error - Compute discrete analogue of the error on the grid.

	properties
		xmin	% Leftmost x-value in rectangular domain.
		xmax	% Rightmost x-value in rectangular domain.
		ymin	% Lowest y-value in rectangular domain.
		ymax	% Highest y-value in rectangular domain.
		nx		% Number of gridpoints in x-direction.
		ny		% Number of gridpoints in y-direction.
		N		% Total number of interior gridpoints.
		hx		% Grid spacing in the x-direction.
		hy		% Grid spacing in the y-direction.
		x		% Vector of x-values.
		y		% Vector of y-values.
		X		% X in [X,Y] = meshgrid(x,y).
		Y		% Y in [X,Y] = meshgrid(x,y)
	end
	
	methods
		
		function self = domain(spatial_data,size_data)
			% Get information from spatial_data
			self.xmin = spatial_data(1);
			self.xmax = spatial_data(2);
			self.ymin = spatial_data(3);
			self.ymax = spatial_data(4);
			
			% Get information from size_data
			self.nx = size_data(1);
			self.ny = size_data(2);
			
			% Calculate quantities which depend on the above
			self.N = self.nx*self.ny;
			self.hx = (self.xmax-self.xmin)/(self.nx+1);
			self.hy = (self.ymax-self.ymin)/(self.ny+1);
			self.x = (self.xmin+self.hx : self.hx : self.xmax-self.hx)';
			assert(length(self.x) == self.nx)
			self.y = (self.ymin+self.hy : self.hy : self.ymax-self.hy)';
			assert(length(self.y) == self.ny)
			[self.X,self.Y] = meshgrid(self.x,self.y);
		end
		
		function M=vec2mat(this,m)
		%VEC2MAT   Reshape a vector to a matrix in accordance with domain.
			M = flipud(reshape(m,this.nx,this.ny)');
		end
		
		function m=mat2vec(this,M)
		%MAT2VEC   Reshape a matrix to a vector in accordance with domain.
			m = reshape(flipud(M)',this.N,1);
		end
		
		function subs=loc2sub(this,locs)
		%LOC2SUB   Transform a location in the domain to its corresponding
		%   subscripts on the grid.
		subs = [round((locs(1,:)-this.xmin)/this.hx);
				round((locs(2,:)-this.ymin)/this.hy)];
		end
		
		function inds = sub2ind(this,subs)
		%SUB2IND   Transform a pair of subscripts on the grid to their
		%   corresponding flat-index in a vector.
			inds = sub2ind([this.nx this.ny],subs(1,:),subs(2,:));
		end
		
		function inds = loc2ind(this,locs)
		%LOC2IND   Transform a location in the domain to its corresponding
		%   flat-index in a vector.
			inds = sort(this.sub2ind(this.loc2sub(locs)));
			if ~(unique(inds) == inds)
				warn('DOMAIN: Non uniqueness of indices!')
				inds = unique_inds;
			end
		end
		
		function b = pt_src(this,xloc,yloc)
		%LOC2IND   Create a point source (discrete Dirac delta) at a given
		%   point  in the domain.
			b = zeros(this.N,1);
			bsub = this.loc2sub([xloc;yloc]);
			bind = this.sub2ind(bsub);
			b(bind) = 1/(this.hx*this.hy);
		end
		
		function [inds,W,IW] = window(this,info)
		%WINDOW   Obtain relevant entities for working with a window in
		%   the domain.
		%
		%	[INDS,W,IW] = WINDOW(THIS,INFO) returns flat indices of the
		%	variables which lie inside of the window, a matrix the size of
		%	the domain which indicates points that are inside of the
		%	window, and a matrix that prolongates variables from the inside
		%	of the window to the rest of the domain.
		%
		%   A window is a region of the domain where unknown values of the
		%   model lie. It is described by the struct INFO, which requires a
		%   field 'type' that either the string disk, rectangle_inner,
		%   rectangle_outer, or all.
		%   
		%   If type is set to disk, fields center and radius of info
		%   describe this disk, and the window is the interior. If type is
		%   set to rectangle_inner or rectangle_outer, the field bounds
		%   should be a vector of the form [xmin xmax ymin ymax] which
		%   describe this rectangle, and the window is the inside or
		%   outside of this rectangle. If type is set to all, the window is
		%   the whole domain (i.e., there is no window).
		
			all_inds = 1:this.N;
			vec = @(X) reshape(X,[],1);
			x = vec(this.X');
			y = vec(this.Y');
			switch info.type
				case 'disk'
					c = info.center;
					r = info.radius;
					inds = all_inds((x-c(1)).^2+(y-c(2)).^2 < r^2);
					W = double((this.X-c(1)).^2+(this.Y-c(2)).^2 < r^2);
				case {'rectangle_inner','rectangle_outer'}
					xm = info.bounds(1);
					xM = info.bounds(2);
					ym = info.bounds(3);
					yM = info.bounds(4);
					if strcmp(info.type,'rectangle_outer')
						inds = all_inds(x <= xm | x >= xM | y <= ym | y >= yM);
						W = double(this.X <= xm | this.X >= xM | this.Y <= ym | this.Y >= yM);
					else
						inds = all_inds(x >= xm & x <= xM & y >= ym & y <= yM);
						W = double(this.X >= xm & this.X <= xM & this.Y >= ym & this.Y <= yM);
					end
				case 'pml'
					wpml = info.width;
					xm = this.xmin+wpml;
					xM = this.xmax-wpml;
					ym = this.ymin+wpml;
					yM = this.ymax-wpml;
					inds = all_inds(x <= xm | x >= xM | y <= ym | y >= yM);
					W = double(this.X <= xm | this.X >= xM | this.Y <= ym | this.Y >= yM);
				case 'all'
					inds = all_inds;
					W = ones(this.ny,this.nx);
			end
			I = speye(this.N);
			IW = I(:,inds);
		end
		
		function imagesc(this,m,meshgrid)
			if ~isreal(m)
				warning('DOMAIN: Discarding complex parts of m.')
				m = real(m);
			end
			if size(m,2) == 1
				if nargin == 2 || meshgrid == false
					imagesc(this.vec2mat(m))
				else
					imagesc(flipud(m))
				end
			else
				if nargin == 2
					imagesc(m);
				end
			end
		end

		function plot(this,x,y,marker,varargin)
			subs = this.loc2sub([x;y]);
			plot(subs(1,:),this.ny-subs(2,:),marker,varargin{:})
		end
		
		function err = error(this,A,B)
			if size(A,2)==1
				err = this.hx*this.hy*norm(A-B);
			else
				err = this.hx*this.hy*norm(A-B,'fro');
			end
		end
		
	end
end