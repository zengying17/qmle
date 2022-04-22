function [g0,grsqall,H0,Hall,phis,chis,gs]=ghfun(lambda,vare,vara,mn,en,mr,er,Xstar,Xbar,Ybar,Ystar,Dvare,DvareR,groupid,blockid)
    varen=Dvare*vare; %n*1 vector of vare.
    varer=DvareR*vare; %r*1 vector of vare 
    Je=size(vare,1);
eJ=ones(Je,1);
kx=size(Xbar,2);
k0=zeros(kx,1);
ki=eye(kx); %kx dimentional identity matrix.

P_o=en./(Dvare*vare+mn*vara);
p_s=((lambda-1)*en+mn)./(mn-en);
beta=(Xstar'*(en./varen.*Xstar)+Xbar'*(P_o.*Xbar))\(Xstar'*(en./varen.*(p_s.*Ystar))+(1-lambda)*Xbar'*(P_o.*Ybar));
ustar=p_s.*Ystar-Xstar*beta;
ubar=Ybar*(1-lambda)-Xbar*beta;
pr_s1=(mr-er)./(mr-er+lambda*er);
qr_s1=1/(1-lambda)*er;

pn_o1=en./varen;
qn_o1=en./(varen+mn*vara);
pr_o1=er./varer;
qr_o1=er./(varer+mr*vara);

pn_w=-en./(mn-en);
qn_w=en;
pr_w=-er./(mr-er);
qr_w=er;

g1=-1*((mr-er)'*(pr_s1.*pr_w)+er'*(qr_s1.*qr_w))+(pn_o1.*pn_w)'*(ustar.*Ystar)+(qn_o1.*qn_w)'*(ubar.*Ybar);
g2=(-0.5*(er'*(qr_o1.*mr))+0.5*((qn_o1.^2.*mn)'*(ubar.^2)));
g3=(-0.5*DvareR'*((mr-er).*pr_o1+qr_o1)+0.5*Dvare'*((pn_o1.^2).*(ustar.^2)+(qn_o1.^2).*(ubar.^2)));

R=max(groupid);
phis=cell(1,R);
chis=cell(1,R);
gs=cell(1,R);
grsqall=zeros(2+Je+kx,2+Je+kx);
for i=1:max(groupid)
    m=mean(mn(groupid==i));
    j=mean(blockid(groupid==i));
    varej=mean(varen(groupid==i));
    ustarr=ustar(groupid==i);
    ubarr=ubar(groupid==i);
    Xstarr=Xstar(groupid==i,:);
    Xbarr=Xbar(groupid==i,:);
phir=[	1/((m-1+lambda)*varej),-m/((1-lambda)*(varej+m*vara)),1/((m-1+lambda)*varej)*beta',-m/((1-lambda)*(varej+m*vara))*beta'; ...
	   	0,-m^2/(2*(varej+m*vara)^2),k0',k0';...
		full(1:Je==j)'.*[-eJ./(2*vare.^2),-m*eJ./(2*(vare.^2+m*vara*eJ).^2),zeros(Je,kx),zeros(Je,kx)];...
	 	k0,k0,-1/varej*ki,-m/(varej+m*vara)*ki]; 
chir=[ustarr'*ustarr-(m-1)*varej;ubarr'*ubarr/m-(vara+varej/m);Xstarr'*ustarr;Xbarr'*ubarr/m];
gr=phir*chir;
grsq=gr*gr';
phis{1,i}=phir;
chis{1,i}=chir;
gs{1,i}=gr;
grsqall=grsqall+grsq;
end
grsqall=(grsqall+grsqall')/2;

h11=-1*((mr-er)'*(pr_s1.^2.*pr_w.^2)+er'*(qr_s1.^2.*qr_w.^2))-(pn_w.^2.*pn_o1)'*(Ystar.^2)-(qn_w.^2.*qn_o1)'*(Ybar.^2);
h12=(-1*(ubar.*Ybar)'*(qn_o1.^2.*qn_w.*mn));
h22=(0.5*(er'*(qr_o1.^2.*mr.^2))-(ubar.*ubar)'*(qn_o1.^3.*mn.^2));

h31=-1*Dvare'*((ustar.*Ystar).*(pn_o1.^2.*pn_w)+(ubar.*Ybar).*(qn_o1.^2.*qn_w));
h33=diag(0.5*DvareR'*((mr-er).*(pr_o1.^2)+(qr_o1.^2))-Dvare'*((ustar.*ustar).*(pn_o1.^3)+(ubar.*ubar).*(qn_o1.^3)));
h32=(0.5*DvareR'*(qr_o1.^2.*mr)-Dvare'*((ubar.*ubar).*(qn_o1.^3.*mn)));

h41=-1*Xstar'*(pn_o1.*pn_w.*Ystar)-Xbar'*(qn_o1.*qn_w.*Ybar);
h42=(-1*Xbar'*(ubar.*qn_o1.^2.*mn));
h34=-1*Dvare'*(Xstar.*(ustar.*(pn_o1.^2))+Xbar.*(ubar.*qn_o1.^2));
h44=-1*Xstar'*(pn_o1.*Xstar)-1*Xbar'*(qn_o1.*Xbar);

H11=[h11,h12,h31'; ...
    h12,h22,h32'; ...
    h31,h32,h33];
H12=[h41';h42';h34];
H21=H12';
H22=(h44+h44')/2;
Hall=[H11,H12;H21,H22];
Hall=(Hall+Hall')/2;
g0=-[g1;g2;g3];
H0=H11-H12/H22*H21;
end