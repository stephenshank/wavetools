function locs=sources_and_receivers(num,info)

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