function [x,out]=lbfgs(CostFunction,x0,opts)
tic;

% implementation of limited memory BFGS algorithm as described in the Byrd
% Nocedal Schnabel 1994 paper (NOT the 1989) "Representation of Quasi-Newton 
% matrices..."
% CostFunction is a function handle that takes one vector argument of the size 
% of x0 and returns 2 variables the (scalar) value of the cost and its gradient.
% x0 is the starting point of the optimization
% opts is a structure with entries defined in the initialization procedure
% below


% Initialization (checks the members of the opts structure ; if not defined
% assign a default value 
if isfield(opts,'Niter') % number of iterations of the optimization
    Niter=opts.Niter;
else
    Niter=100;
end
if isfield(opts,'optsLineSearch')% structure containing the options for the linesearch (see LineSearchWolfe.m)
    optsLineSearch=opts.optsLineSearch;
else
    optsLineSearch={};
end
% if isfield(opts,'init_t')% initial guess for the step size in the line search
%     init_t=opts.init_t;
% else
%     init_t='';
% end
if isfield(opts,'mem')% memory (number of former steps used to build the approximate inverse Hessian)
    mem=opts.mem;
else
    mem=4;
end
if isfield(opts,'nCallMax')% maximum allowed number of calls of CostFunction (stop when reached)
    nCallMax=opts.nCallMax;
else
    nCallMax=Inf;
end
if isfield(opts,'TMax')% maximum allowed time (stop when reached)
    TMax=opts.TMax;
else
    TMax=Inf;
end
if isfield(opts,'tol')% tolerance over the reduction of the norm of the gradient
    tol=opts.tol;
else
    tol=1e-5;
end



J=zeros(Niter,1);% stores the values of the cost
ng=zeros(Niter,1);% stores the values of the norm of the gradient 
callPerIter=zeros(Niter,1);% stores the number of inner calls (line search) for each outer iteration
hist=zeros(size(x0,1),Niter);% stores the optimization variable for each step
step=zeros(Niter,1);% stores the step sizes
time=zeros(Niter,1);% store the time for each iteration



% intialize l-BFGS
n=size(x0,1);
x=x0;
S=zeros(n,0);
Y=zeros(n,0);


% the initial guess for H (inverse Hessian is gamma*I with gamma=J/||gradJ||^2, linear approximation estimation)
[Jx,gFx]=CostFunction(x);
J0=Jx;
ng0=norm(gFx);
gamma=Jx/ng0^2;


nCall=1;



for i=1:Niter
    t1=toc;
    
    % builds the approximate inverse Hessian

    W=[S gamma*Y];    
    c1=W'*(-gFx);
    R=S'*Y;
    R=tril(R);
    D=diag(diag(R));
    iR=inv(R);
    M=[iR'*(D+gamma*(Y'*Y))*iR -iR';-iR zeros(size(R))];
    c2=M*c1;
    c3=W*c2;
    
    % get the direction (H*-gradJ)
    p=-gamma*gFx+c3;
    
    % linesearch in the direction p
    [t,xp,Jxp,gFxp,nInnerCall]=LineSearchVincent(x,p,Jx,gFx,CostFunction,optsLineSearch);
    step(i)=t;
    if nInnerCall>7
        fprintf('Lots of line search loops\n');
    end
    
    callPerIter(i)=nInnerCall;  
    nCall=nCall+nInnerCall;
    % update the l-BFGS objects
    s=xp-x;
    y=gFxp-gFx;
    
   
    if i<=mem
        S=[S s];
        Y=[Y y];
    else
        S=circshift(S,[0 -1]);
        S(:,mem)=s;
        Y=circshift(Y,[0 -1]);
        Y(:,mem)=y;
    end

    gamma=(y'*s)/(y'*y);

      
    ng(i)=norm(gFxp);
    
    xlast=x;
    plast=p;
    
    x=xp;
%     save('profile','x')
    Jx=Jxp;
    gFx=gFxp;

    hist(:,i)=x;
    J(i)=Jxp;
%     save('misfit','J')
    
    
    
    t2=toc;
    time(i)=t2;
    fprintf('Iteration %d, Jk/J0=%.3d time %.2fs ',i,J(i)/J0,t2-t1);
    fprintf('|gk|/|g0| = %.2d\n',ng(i)/ng0);
    if ng(i)<tol*ng0
        fprintf('small gradient reached |g|<%.1d|g0|\n',tol);
        break;
    end
    if nCall>nCallMax
        fprintf('maximum number of calls of cost function /gradient reached : nCall=%d\n',nCall);
        break;
    end
    if t2>TMax
        fprintf('maximum time reached : T=%.2fs\n',t2);
        break;
    end
    


    
    
    
end
if i==Niter
    fprintf('Maximum number of iterations reached Niter=%d\n',Niter);
end    

% returns the optimization details in the structure out
out.time=time;
out.J=J;
out.callPerIter=callPerIter;
out.ng=ng;
out.xlast=xlast;
out.plast=plast;
out.hist=hist;
out.nCall=nCall;
out.step=step;

function [test,xt,qt,gt,nInnerCall]=LineSearchVincent(x,d,q0,g0,CostFunction,opts)    

% implements linesearch with Wolfe's conditions 
% x is the point at which the line-search is required
% d is the direction of the linesearch
% q0 is the value of the function at x
% g0 is the gradient of the function at x
% CostFunction is a function handle that takes one vector argument of the size 
% of x0 and returns 2 variables the (scalar) value of the cost and its gradient.
% opts is a structure (see initialization)

% Initialization
if isfield(opts,'m1')% parmater for 1st Wolfe's condition
    m1=opts.m1;
else
    m1=0.001;
end
if isfield(opts,'m2')% parmater for 2nd Wolfe's condition
    m2=opts.m2;
else
    m2=0.9;
end
if isfield(opts,'t0')% initial guess for the step
    t0=opts.t0;
else
    t0=1;
end
if isfield(opts,'interp')% type of interpolation ('simple' or ' cubic') ' cubic' has not been tested in a long time, I would advise to stick with simple
    interp=opts.interp;
else
    interp='simple';
end
if isfield(opts,'extrap')% type of extrapolation ('simple' or ' cubic') ' cubic' has not been tested in a long time, I would advise to stick with simple
    extrap=opts.extrap;
else
    extrap='simple';
end
if isfield(opts,'maxt')% maximum step allowed
   maxt=opts.maxt;
else
    maxt=1e12;
end


a=2;
theta=0.1;
% mint=eps;

test=t0;
tolStep=1e-2;

tL=0;tR=0;
qprime0=g0'*d;
qt=Inf;
qprimet=0;
comptLineSearch=0;
nInnerCall=0;
tm=0;
qtm=q0;
t=t0;
qprimetm=qprime0;

go=((qt>q0+m1*t*qprime0 || qprimet<m2*qprime0));
while go
    test=t;
    xt=x+t*d;
    [qt,gt]=CostFunction(xt);% computes cost and gradient at x+t*d
    nInnerCall=nInnerCall+1;
    qprimet=gt'*d;
    if qt>q0+m1*t*qprime0
        tR=t;
    else
        tL=t;
    end
    
    if tR==0
        if strcmp(extrap,'simple')
            t=a*t;
        elseif strcmp(extrap,'cubic')
            
            tp=CubicFitting(t,tm,qt,qtm,qprimet,qprimetm);
            tm=t;
            qtm=qt;
            qprimetm=qprimet;
            t=max(tp,a*t);
        else
            fprintf('error : Unknown extrapolation option\n');
            break;
        end
        
    else
        if strcmp(interp,'simple')
            t=(tR+tL)/2;
        elseif strcmp(interp,'cubic')
            tp=CubicFitting(t,tm,qt,qtm,qprimet,qprimetm);
            tm=t;
            qtm=qt;
            qprimetm=qprimet;
            if mod(comptLineSearch,2)==0
                t=min(tp,tR-theta*(tR-tL));
            else
                t=max(tp,tL+theta*(tR-tL));
            end

        else
            fprintf('error : Unknown extrapolation option\n');
            break;
        end

    end
    comptLineSearch=comptLineSearch+1;
    
    if ~(qt>q0+m1*t*qprime0 || qprimet<m2*qprime0)
        fprintf('Found a step that satisfies Wolfe''s conditions\n');
        go=0;
    elseif tL>maxt
        fprintf('Forced stop, tL is too large\n');
        go=0;
    elseif (abs((tR-tL)/max(tR,tL))<tolStep)
        fprintf('Forced stop, tL and tR are too close\n');
        go=0;

    elseif comptLineSearch>=30
        fprintf('Forced stop, too many inner loops\n');
        go=0;
    end
        
    if go
        fprintf('Inner loop %d tL=%d tR=%d t=%d\n',comptLineSearch,tL,tR,t);
    end
end
if comptLineSearch==0
    fprintf('NO LINE SEARCH?? q''(0)=%d\n',qprime0);
end
fprintf('Linesearch over : %d evaluation(s), estimated step=%.2d\n',comptLineSearch,test);