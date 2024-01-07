%% steatdy-soil water-solution
function [theta_s]=steadysolution(Nlay,omega,K_s,n_0,alpha,q_s,q_p,L)
k_=zeros(1,Nlay);
Kc=zeros(1,Nlay);
theta_s=zeros(1,L(Nlay)/0.01);
theta_p=zeros(1,L(Nlay)/0.01);
K_c=zeros(1,L(Nlay)/0.01);
K_p=zeros(1,L(Nlay)/0.01);
qqp=zeros(1,L(Nlay)/0.01);
lamber1=zeros(1,Nlay);
lamber2=zeros(1,Nlay);
for i_lay=1:Nlay
k_(i_lay)=K_s(i_lay)/n_0(i_lay); %Parameter
lamber1(i_lay)=0.5-0.5.*(1+(4.*1i.*omega)./(k_(i_lay).*alpha(i_lay))).^(0.5); %Constant  
lamber2(i_lay)=0.5+0.5.*(1+(4.*1i.*omega)./(k_(i_lay).*alpha(i_lay))).^(0.5); %Constant  
end

Kirchhoff1=zeros(1,Nlay);
Kirchhoff2=zeros(1,Nlay);
fluxcoef1=zeros(1,Nlay);
fluxcoef2=zeros(1,Nlay);
fKirchhoff1=zeros(1,Nlay-1);
fKirchhoff2=zeros(1,Nlay-1);
ffluxcoef1=zeros(1,Nlay-1);
ffluxcoef2=zeros(1,Nlay-1);

Kc1=zeros(1,Nlay);
Kc2=zeros(1,Nlay-1);


for i_lay=1:Nlay
%     Kirchhoff1(i_lay)=exp(lamber1(i_lay)*alpha(i_lay)*L(i_lay));
%     Kirchhoff2(i_lay)=exp(lamber2(i_lay)*alpha(i_lay)*L(i_lay));
%     fluxcoef1(i_lay)=(1-lamber1(i_lay))*Kirchhoff1(i_lay);
%     fluxcoef2(i_lay)=(1-lamber2(i_lay))*Kirchhoff2(i_lay);
    Kc1(i_lay)=exp(alpha(i_lay)*L(i_lay))/alpha(i_lay);
end
for i_lay=1:(Nlay-1)
%     fKirchhoff1(i_lay)=exp(lamber1(i_lay+1)*alpha(i_lay+1)*L(i_lay));
%     fKirchhoff2(i_lay)=exp(lamber2(i_lay+1)*alpha(i_lay+1)*L(i_lay));
%     ffluxcoef1(i_lay)=(1-lamber1(i_lay+1))*fKirchhoff1(i_lay);
%     ffluxcoef2(i_lay)=(1-lamber2(i_lay+1))*fKirchhoff2(i_lay);
    Kc2(i_lay)=exp(alpha(i_lay+1)*L(i_lay))/alpha(i_lay+1);
end
BC=q_p;
DE=q_s;
% AA=zeros(2*Nlay,2*Nlay);
%BB=zeros(2*Nlay,1);
% AA(1,1)=(1-lamber1(1));
% AA(1,1+Nlay)=(1-lamber2(1));
DD=zeros(2*Nlay,2*Nlay);
DD(1,1+Nlay)=1;
EE=zeros(2*Nlay,1);
for i=1:Nlay-1
%     AA(i+1,i)=fluxcoef1(i);
%     AA(i+1,i+Nlay)=fluxcoef2(i);
%     AA(i+1,i+1)=-ffluxcoef1(i);
%     AA(i+1,i+1+Nlay)=-ffluxcoef2(i);
%     AA(i+Nlay,i)=Kirchhoff1(i)/K_s(i);
%     AA(i+Nlay,i+Nlay)=Kirchhoff2(i)/K_s(i);
%     AA(i+Nlay,i+1)=-fKirchhoff1(i)/K_s(i+1);
%     AA(i+Nlay,i+1+Nlay)=-fKirchhoff2(i)/K_s(i+1);
   % DD(i+1,i)=Kc1(i)/K_s(i);
    %DD(i+1,i+Nlay)=-Kc2(i)/K_s(i+1);
    %DD(i+1,i+1)=1/K_s(i);
    %DD(i+1,i+1+Nlay)=-1/K_s(i+1);
   % DD(i+Nlay,i+1)=1;
   % DD(i+Nlay,i+1+Nlay)=-1;
   DD(i+1,i+Nlay)=1;
   DD(i+1,i+1+Nlay)=-1;
   DD(i+Nlay,i)=Kc1(i)/K_s(i);
   DD(i+Nlay,i+Nlay)=1/K_s(i);
   DD(i+Nlay,i+1)=-Kc2(i)/K_s(i+1);
   DD(i+Nlay,i+1+Nlay)=-1/K_s(i+1);
end
 %AA(2*Nlay,Nlay)=Kirchhoff1(Nlay);
%AA(2*Nlay,2*Nlay)=Kirchhoff2(Nlay);
DD(2*Nlay,Nlay)=1;
EE(1)=DE;
%BB(1)=BC;
%CC=AA\BB;
GG=DD\EE;

deep=linspace(0.00,L(Nlay),L(Nlay)/0.01+1);
for kk=1:(L(1)/0.01+1)
    K_c(kk)=GG(1)/alpha(1)*exp(alpha(1)*deep(kk))+GG(1+Nlay);
%    K_p(kk)=CC(1)*exp(lamber1(1)*alpha(1)*deep(kk))+CC(1+Nlay)*exp(lamber2(1)*alpha(1)*deep(kk));
%    qqp(kk)=(1-lamber1(1))*CC(1)*exp(lamber1(1)*alpha(1)*deep(kk))+(1-lamber2(1))*CC(1+Nlay)*exp(lamber2(1)*alpha(1)*deep(kk));
    theta_s(kk)=n_0(1).*(K_c(kk)/K_s(1));
%    theta_p(kk)=n_0(1).*(K_p(kk)/K_s(1));
end
    for i_lay=2:Nlay
        for kk=1:(L(Nlay)/0.01+1)
        if deep(kk)<=L(i_lay)&&deep(kk)>L(i_lay-1)
        K_c(kk)=GG(i_lay)/alpha(i_lay)*exp(alpha(i_lay)*deep(kk))+GG(i_lay+Nlay);
%        K_p(kk)=CC(i_lay)*exp(lamber1(i_lay)*alpha(i_lay)*deep(kk))+CC(i_lay+Nlay)*exp(lamber2(i_lay)*alpha(i_lay)*deep(kk));
%        qqp(kk)=(1-lamber1(i_lay))*CC(i_lay)*exp(lamber1(i_lay)*alpha(i_lay)*deep(kk))+(1-lamber2(i_lay))*CC(i_lay+Nlay)*exp(lamber2(i_lay)*alpha(i_lay)*deep(kk));
        theta_s(kk)=n_0(i_lay).*(K_c(kk)/K_s(i_lay));
%        theta_p(kk)=n_0(i_lay).*(K_p(kk)/K_s(i_lay));
        end
        end
    end
