function loc=sources_2d(num,info)

loc=zeros(num,2);
switch info.type
	case 'circle'
		theta=2*pi*(0:num-1)/num;
		loc(:,1)=info.center(1)+info.radius*cos(theta);
		loc(:,2)=info.center(2)+info.radius*sin(theta);
	case 'hline'
		loc(:,1)=linspace(info.left,info.right,num);
		loc(:,2)=info.height;
	case 'vline'
		loc(:,1)=info.pos;
		loc(:,2)=linspace(info.top,info.bottom,num);
end