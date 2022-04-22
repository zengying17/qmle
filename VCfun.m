function [Ups,Gam]=VCfun(Xbar,Xstar,N,mn,mr,lambda,vare,vara,beta,mu3e,mu3a,mu4e,mu4a,Rtab,blockid)
Je=size(vare,1);
eJ=ones(Je,1);
kx=size(Xbar,2);
k0=zeros(kx,1);
ki=eye(kx); %kx dimentional identity matrix.
Ups=zeros(2+Je+kx);
Gam=zeros(2+Je+kx);
for j=1:Je
for m=2:max(mr)
Rmj=sum(Rtab.mn==m&Rtab.blockid==j); %number of groups with size equals to m
if Rmj>0
Dm=Xstar(mn==m&blockid==j,:);
Dm2=Dm'*Dm;
Bm=Xbar(mn==m&blockid==j,:);
Bm2=Bm'*Bm;
Cm=sum(Bm,1)/m;

varej=vare(j);
mu3ej=mu3e(j);
mu4ej=mu4e(j);
phimj=[	1/((m-1+lambda)*varej),-m/((1-lambda)*(varej+m*vara)),1/((m-1+lambda)*varej)*beta',-m/((1-lambda)*(varej+m*vara))*beta'; ...
	   	0,-m^2/(2*(varej+m*vara)^2),k0',k0';...
		full(1:Je==j)'.*[-eJ./(2*vare.^2),-m*eJ./(2*(vare+m*vara*eJ).^2),zeros(Je,kx),zeros(Je,kx)];...
	 	k0,k0,-1/varej*ki,-m/(varej+m*vara)*ki]; 
	 	 
PsiGmj=blkdiag(2*(m-1)*varej^2*Rmj/N,2*(vara+varej/m)^2*Rmj/N,varej*Dm2/N,(vara+varej/m)*(Bm2/m)/N);

Psi11=Rmj/N*blkdiag(2*(m-1)*varej^2,2*(vara+varej/m)^2+(mu4a-3*vara^2))+Rmj/N*(mu4ej-3*varej^2)*[(m-1)^2/m,(m-1)/m^2;(m-1)/m^2,1/m^3];
Psi21=[k0,k0;(m-1)/m*mu3ej*Cm'/N,(mu3a+mu3ej/m^2)*Cm'/N];
Psi22=blkdiag(varej*Dm2/N,(vara+varej/m)*(Bm2/m)/N);
Psimj=[Psi11,Psi21';Psi21,Psi22];
Upsmj=phimj*Psimj*phimj';
Gammj=phimj*PsiGmj*phimj';

Ups=Ups+Upsmj;
Gam=Gam+Gammj;
end
end
Ups=(Ups+Ups')/2;
Gam=(Gam+Gam')/2;
end