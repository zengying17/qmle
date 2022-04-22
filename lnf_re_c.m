function [f,g,H]=lnf_re_c(x,N,R,mn,mr,en,er,Ystar,Ybar,Xstar,Xbar,Dvare,DvareR)
    lambda=(1-exp(x(1)))/(1+exp(x(1)));
    vara=x(2)^2;
    vare=(x(3:end).^2)'; %Je*1 vector of vare.
    varen=Dvare*vare; %n*1 vector of vare.
    varer=DvareR*vare; %r*1 vector of vare 

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
H21=[h41,h42,h34'];
H22=h44;

f=0.5*N*log(2*pi())-R*log(abs(1-lambda))-(mr-er)'*log(er+lambda*er./(mr-er))  ...
+0.5*mr'*log(varer)+0.5*er'*log(er+vara*(mr./varer)) ...
+0.5*(ustar'*(1./varen.*ustar)+ubar'*(P_o.*ubar));  

der1=-2*exp(x(1))/(1+exp(x(1)))^2;
der2=2*exp(x(1))*(exp(x)-1)*(exp(x(1))+1)/(1+exp(x(1)))^4;
P=diag([der1,2*x(2),2*x(3:end)]);
if nargout >1
g0=-[g1;g2;g3];
g=P*g0;
if nargout >2
H=-(P*H11*P-P*H12/H22*H21*P)+diag([der2,2*ones(1,size(x,2)-1)]'.*g0); 
end
end