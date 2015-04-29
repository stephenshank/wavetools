classdef domain
	
	properties
		xm
		xM
		ym
		yM
		nx
		ny
		N
		hx
		hy
		x
		y
		X
		Y
	end
	
	methods
		
		function self = domain(spatial_data, size_data)
			self.xm = spatial_data(1);
			self.xM = spatial_data(2);
			self.ym = spatial_data(3);
			self.yM = spatial_data(4);
			self.nx = size_data(1);
			self.ny = size_data(2);
			self.N = self.nx*self.ny;
			self.hx = (self.xM-self.xm)/(self.nx+1);
			self.hy = (self.yM-self.ym)/(self.ny+1);
			self.x = (self.xm+self.hx : self.hx : self.xM-self.hx)';
			assert(length(self.x) == self.nx)
			self.y = (self.ym+self.hy : self.hy : self.yM-self.hy)';
			assert(length(self.y) == self.ny)
			[self.X,self.Y] = meshgrid(self.x,self.y);
		end
		
		function M=m2M(this,m)
			M = flipud(reshape(m,this.nx,this.ny)');
		end
		
		function m=M2m(this,M)
			m = reshape(flipud(M)',this.N,1);
		end
		
		function subs=loc2sub(this,locs)
			subs = [round((locs(1,:)-this.xm)/this.hx);
					round((locs(2,:)-this.ym)/this.hy)];
		end
		
		function inds = sub2ind(this,subs)
			inds = sub2ind([this.nx this.ny],subs(1,:),subs(2,:));
		end
		
		function inds = loc2ind(this,locs)
			inds = sort(this.sub2ind(this.loc2sub(locs)));
			if ~(unique(inds) == inds)
				warn('DOMAIN: Non uniqueness of indices!')
				inds = unique_inds;
			end
		end
		
		function b = pt_src(this,xloc,yloc)
			b = zeros(this.N,1);
			bsub = this.loc2sub([xloc;yloc]);
			bind = this.sub2ind(bsub);
			b(bind) = 1/(this.hx*this.hy);
		end
		
		function [inds,W,IW] = window(this,info)
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
				case 'rectangle'
					xm = info.bounds(1);
					xM = info.bounds(2);
					ym = info.bounds(3);
					yM = info.bounds(4);
					inds = all_inds(x <= xm | x >= xM | y <= ym | y >= yM);
					W = double(this.X <= xm | this.X >= xM | this.Y <= ym | this.Y >= yM);
				case 'pml'
					wpml = info.width;
					xm = this.xm+wpml;
					xM = this.xM-wpml;
					ym = this.ym+wpml;
					yM = this.yM-wpml;
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
				warning('MESH: Discarding complex parts of m.')
				m = real(m);
			end
			if size(m,2) == 1
				if nargin == 2 || meshgrid == false
					imagesc(this.m2M(m))
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
			err = this.hx*this.hy*norm(A-B,'fro');
		end
		
	end
end