clear all
close all
clc

f = 50;			% frequency
n = 10*f;		% number of grid points
omega = 2*pi*f;	% angular frequency
x = .2;			% x coordinate of point source
y = .8;			% x coordinate of point source

% Free space solution
A = helmholtz_2d(ones(n),f,n);
b = pt_src(n,x,y);
U = reshape(A\b,n,n);

figure
subplot(3,2,1)
[X,Y]=meshgrid(linspace(0,1,n));
Utrue=flipud(1i*besselh(0,omega*abs(X+1i*Y-(x+1i*y)))/4);
imagesc(real(Utrue));

subplot(3,2,2)
imagesc(flipud(real(U)'))
cvec=caxis;

% Shepp-Logan
ctr = 2;
c = 1 + ctr*phantom(n);
A = helmholtz_2d(1./c.^2,f,n);
b = pt_src(n,x,y);
U = reshape(real(A\b),n,n);

subplot(3,2,3)
imagesc(c)
hold on
plot(round((n+1)*x),n-round((n+1)*y),'r.','MarkerSize',10)
hold off

subplot(3,2,4)
imagesc(flipud(U'))

% Marmousi
c = marmousi(n);
A = helmholtz_2d(1./c.^2,f,n);
U = reshape(real(A\b),n,n);

subplot(3,2,5)
imagesc(c)
hold on
plot(round((n+1)*x),n-round((n+1)*y),'r.','MarkerSize',10)
hold off

subplot(3,2,6)
imagesc(flipud(U'))

titles={'Free space solution (analytic)','Computed free space solution', ...
	'Shepp-Logan (medical)','Computed Shepp-Logan solution', ...
	'Marmousi (geophysics)','Computed Marmousi solution',};
for i=1:6
	subplot(3,2,i)
	if ismember(i,[1 2 4 6])
		caxis(cvec)
	end
	set(gca, 'XTick', []);
	set(gca, 'YTick', []);
	title(titles{i},'fontweight','bold','fontsize',14)
end
set(gcf,'Position',[50 120 480 660])