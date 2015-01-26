function A = helmholtz_fos_1d(m,f)

n=length(m)-2;
omega=2*pi*f;
h=1/(n+1);
e=ones(n+2,1);
A=spdiags([-e 2*e -e],-1:1,n+2,n+2)/h^2-omega^2*spdiags(m,0,n+2,n+2);
A(1,1)=(1/h-1i*omega);
A(1,2)=-1/h;
A(n+2,n+1)=-1/h;
A(n+2,n+2)=(1/h-1i*omega);