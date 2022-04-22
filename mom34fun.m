%This is the function for estimating the third and fourth moments
%[mu3e,mu4e,mu3a,mu4a]=mom34fun(ubar,ustar,mn,mr,vare,vara,Dvare);
%2019-02-05, demean the u. 
function [mue3,mue4,mua3,mua4]=mom34fun(ubar,ustar,mn,mr,vare,vara,Dvare,DvareR)
ubar=ubar-mean(ubar); %This is optional
udot=ustar;
N=size(mn,1);
en=ones(N,1);
R=size(mr,1);

ubar3=ubar.^3;
ubar4=ubar.^4;
udot3=udot.^3;
udot4=udot.^4;
udot2bar=udot.^2.*ubar;
fe3a=((udot3./(en-3*en./mn+2*en./(mn.^2)))./mn);
fe3a(isinf(fe3a))=0;
fe3b=((udot2bar./(en./mn-en./(mn.^2)))./mn);
fe3b(isinf(fe3b))=0;

fe3=full(fe3a.*(mn>=3)+fe3b.*(mn<=2));
fe4=((mn.^3)./(mn.^3-4*mn.^2+6*mn-3*en)).*(udot4./mn-3*(mn-en).*(2*mn-3*en)./(mn.^4).*((Dvare*vare).^2));
fe4=full(fe4);
fa3=ubar3./mn-fe3./(mn.^2);
fa4=ubar4./mn-fe4./(mn.^3)-3*(mn-1)./(mn.^4).*(Dvare*vare).^2-6./(mn.^2).*(Dvare*vare)*vara;

mue3=(Dvare'*fe3)./sum(full(DvareR))';
mue4=(Dvare'*fe4)./sum(full(DvareR))';
mua3=mean(fa3)*(N/R);
mua4=mean(fa4)*(N/R);
end