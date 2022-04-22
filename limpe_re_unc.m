function [result,sample]=limpe_re_unc(datatable,yvar,groupvar,xivar,xpvar,hetvar,options,x0)
disp('limpe_re_unc: linear in means peer effects model with random group effects under no constraints')
%% Select the Sample

sampletb=datatable(:,[yvar,groupvar,xivar,xpvar,hetvar]);
TP=ismissing(sampletb);
sampletb=sampletb(~any(TP,2),:); %DTB is the final data table to be used. delete rows with missing value.

tempid=grp2idx(sampletb.(groupvar));
for j=1:max(tempid)
    mnj=sum(tempid==j);
    sampletb.mn(tempid==j)=mnj;
end

T=sampletb(sampletb.mn>1,:); %delete classes with only one student.
T=sortrows(T,groupvar);

%% Define the inputs
y=T.(yvar);
mn=T.mn;
groupid=grp2idx(T.(groupvar));
if isempty(hetvar)==1
    blockid=ones(size(y,1),1);
else
blockid=grp2idx(table2array(T(:,hetvar)));
end
Dvare=sparse(dummyvar(blockid));
Je=size(Dvare,2); % number of categories of vare.
if isempty(xivar)==1
    Xi=[];
else 
    Xi=table2array(T(:,xivar));
end
if isempty(xpvar)==1
    Xp=[];
else 
    Xp=table2array(T(:,xpvar));
end
%% 
if nargin<=6 || isempty(options) ==1
options = optimoptions(@fminunc,'Algorithm','Trust-region','SpecifyObjectiveGradient',true,'DerivativeCheck','off','FinDiffType','Central','Display','final','FunctionTolerance',1e-12,'StepTolerance',1e-12,'MaxFunctionEvaluations',2000,'MaxIterations',2000);
end
if nargin<=7 || isempty(x0) ==1
x0=[0,var(y),var(y)*ones(1,Je)];
end

%% Define matrices
N=size(groupid,1);
R=max(groupid);
en=ones(N,1);
er=ones(R,1);
ximean=zeros(N,size(Xi,2));
xpmean=zeros(N,size(Xp,2));
Ybar=zeros(N,1);

for i=1:max(groupid)
    mi=sum(groupid==i);
    if isempty(xivar)==0 
    ximean(groupid==i,:)=repmat(mean(Xi(groupid==i,:)),mi,1);
    end
    if isempty(xpvar)==0
        xpmean(groupid==i,:)=repmat(mean(Xp(groupid==i,:)),mi,1);
    end
    Ybar(groupid==i,:)=repmat(mean(y(groupid==i,:)),mi,1);
end

Ystar=y-Ybar;
Xstar=zeros(N,1);
Xbar=ones(N,1);
if isempty(Xp)==0      %If we allow for exogeneous peer effects, then X should include the spatial lags of these variables.
    Xbar=[xpmean,Xbar];
    Xstar=[-(en./(mn-en)).*(Xp-xpmean),Xstar];
end
if isempty(Xi)==0
    Xstar=[Xi-ximean,Xstar];
    Xbar=[ximean,Xbar];
end

[~, jb]= rref(Xstar+Xbar);      %index for rows that are not colinear
Xstar= Xstar(:, jb);
Xbar= Xbar(:, jb);
xnames0=[xivar strcat('ave_',xpvar(:))' 'cons'];
xnames=xnames0(:,jb);
estnames=['lambda' 'vara' strcat('vare',cellstr(string(1:Je))) xnames];

Rtab=unique(table(groupid,mn,blockid,Dvare));
mr=Rtab.mn;
DvareR=Rtab.Dvare;

%% Optimization
format long
xinit=x0*NaN;
xinit(1)=log((1-x0(1))/(1+x0(1)));
xinit(2)=x0(2)^0.5;
xinit(3:end)=x0(3:end).^0.5;

% first step: constrained optimization.
[x,lnf,exitflag,out,grad,hessian]=fminunc(@(x) lnf_re_c(x,N,R,mn,mr,en,er,Ystar,Ybar,Xstar,Xbar,Dvare,DvareR),xinit,options);
lambda=(1-exp(x(1)))/(1+exp(x(1)));
vara=x(2)^2;
vare=(x(3:end).^2)'; %Je*1 vector of vare.
xfirst=[lambda,vara,vare'];

[x,lnf,exitflag,out,grad,hessian]=fminunc(@(x) lnf_re_unc(x,N,R,mn,mr,en,er,Ystar,Ybar,Xstar,Xbar,Dvare,DvareR),xfirst,options);
if exitflag<=0
    x0new=x;
    myopt=optimoptions(@fminunc,'Display','final');
[x,lnf,exitflag,out,grad,hessian]=fminunc(@(x) lnf_re_unc(x,N,R,mn,mr,en,er,Ystar,Ybar,Xstar,Xbar,Dvare,DvareR),x0new,myopt);
end

lambda=x(1);
vara=x(2);
vare=(x(3:end))';
varen=Dvare*vare; %n*1 vector of vare.
lnf=-lnf; %Becareful, we were minimizing -lnf.

P_o=en./(varen+mn*vara);
p_s=((lambda-1)*en+mn)./(mn-en);

beta=(Xstar'*((1./varen).*Xstar)+Xbar'*(P_o.*Xbar))\(Xstar'*(1./varen.*p_s.*Ystar)+(1-lambda)*Xbar'*(P_o.*Ybar));
ustar=p_s.*Ystar-Xstar*beta;
ubar=Ybar*(1-lambda)-Xbar*beta;

[mu3e,mu4e,mu3a,mu4a]=mom34fun(ubar,ustar,mn,mr,vare,vara,Dvare,DvareR);
mom34_est=array2table([mu3e',mu4e',mu3a,mu4a],'VariableNames',[strcat('mu3e',cellstr(string(1:Je))),strcat('mu4e',cellstr(string(1:Je))),'mu3a','mu4a']);

[Ups,Gam]=VCfun(Xbar,Xstar,N,mn,mr,lambda,vare,vara,beta,mu3e,mu3a,mu4e,mu4a,Rtab,blockid);
V1=inv(Gam)/N;
V2=N*(V1*Ups*V1);
[g0,grsqall,H0,Hall,phis,chis,gs]=ghfun(lambda,vare,vara,mn,en,mr,er,Xstar,Xbar,Ybar,Ystar,Dvare,DvareR,groupid,blockid);
V3=inv(-Hall);
V4=inv(Hall/grsqall*Hall);

b=[lambda,vara,vare',beta'];
se1=sqrt(diag(V1))';
se2=sqrt(diag(V2))';
se3=sqrt(diag(V3))';
se4=sqrt(diag(V4))';

p1=2*(1-cdf('normal',abs(b./se1),0,1));
p2=2*(1-cdf('normal',abs(b./se2),0,1));
p3=2*(1-cdf('normal',abs(b./se3),0,1));
p4=2*(1-cdf('normal',abs(b./se4),0,1));

resmat=[b;se1;p1;se2;p2;se3;p3;se4;p4];
restab1=array2table(resmat,'VariableNames',estnames,'rownames',{'est','se1','p1','se2','p2','se3','p3','se4','p4'});
restab2=[table(lnf,N,R,exitflag),mom34_est];

bnames=strcat(estnames,'_b');
se1names=strcat(estnames,'_se1');
p1names=strcat(estnames,'_p1');
se2names=strcat(estnames,'_se2');
p2names=strcat(estnames,'_p2');
se3names=strcat(estnames,'_se3');
p3names=strcat(estnames,'_p3');
se4names=strcat(estnames,'_se4');
p4names=strcat(estnames,'_p4');
resvec=[splitvars(table(b,se1,p1,se2,p2,se3,p3,se4,p4,lnf,N,R,exitflag),{'b','se1','p1','se2','p2','se3','p3','se4','p4'},'NewVariableNames',{bnames,se1names,p1names,se2names,p2names,se3names,p3names,se4names,p4names}),mom34_est];
format short
disp(restab1)
disp(restab2)   

%% Export Restults.
result.lambda=lambda;
result.vare=vare;
result.vara=vara;
result.beta=beta;
result.b=b;
result.xnames=xnames;
result.estnames=estnames;
result.restab1=restab1;
result.restab2=restab2;
result.resvec=resvec;
result.lnf=lnf;  
result.N=N;
result.R=R;
result.mom34_est=mom34_est;
result.Ups=Ups;
result.Gam=Gam;
result.grsqall=grsqall;
result.Hall=Hall;
result.V1=V1;
result.V2=V2;
result.V3=V3;
result.V4=V4;
result.x0=x0;
result.xinit=xinit;
result.xfirst=xfirst;
result.x=x;
result.exitflag=exitflag;
result.output=out;
sample.X=Xstar+Xbar;
sample.y=y;
sample.groupid=groupid;
sample.blockid=blockid;
sample.ustar=ustar;
sample.ubar=ubar;
sample.mn=mn;
end