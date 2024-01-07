%% Inversion of assimilation PS
function [ksrenew,ksrmse]=inverse_ps(ks,Nfobs,numbf,Ndep,numdep,Njtime,Lx,K_s,q_s,q_p,alpha,n_0)
load true.mat 
Am0=angm0;
Nlay=20;
L=[0.10,0.20,0.30,0.40,0.50,0.60,0.70,0.80,0.90,1.00,1.10,1.20,1.30,1.40,1.50,1.60,1.70,1.80,1.90,2.00];
num_pam=linspace(1,Nlay,Nlay);
delect_pam=ceil(numdep/10);
lambda=1*ones(1,Njtime);
lamb_adj=1.2;
lgtt=linspace(-3,1,1001);
f=q_switch(0.31,alpha(1),n_0(1),lgtt);
NN=length(f);
for j=1:NN
qp(j)=complex(q_p*cos(0),q_p*sin(0));
end 

%% 

ks=log(ks');
ksrenew=zeros(Nlay,Njtime+1);
ksrenew(:,1)=ks;
refAmaa=zeros(Nfobs,Ndep);
for i=1:(Nfobs)
for    j=1:Ndep
    if(logth(numbf(i),numdep(j)))
refAmaa(i,j)=angm0(numbf(i),numdep(j));
    end
end 
end
refA=refAmaa';
refA=refA(:);
refA(all(refA==0,2),:)=[];
Ndd=size(refA,1);
%% 
delta=0.00001;
for jj=1:Njtime
for i=1:Nfobs
[~,comres(i,:),~]=fresolution(Nlay,2*pi*f(numbf(i)),exp(ksrenew(:,jj)),n_0,alpha,q_s,qp(numbf(i)),L);%ѡȡƵ�ʽ��м���
end
Am1=angle(comres);
virtAma=zeros(Nfobs,Ndep);
for k=1:Ndep 
virtAma(:,k)=Am1(:,numdep(k));
end
virta=zeros(Nfobs,Ndep);
for i=1:(Nfobs)
for    j=1:Ndep
if(logth(numbf(i),numdep(j)))
 virta(i,j)=virtAma(i,j);
end
end
end
virta=virta';
virta=virta(:);
virta(all(virta==0,2),:)=[];
rmse(1)=sqrt(sum((virta-refA).^2)/(Ndd));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
kkk=1;
for kk=1:(Nlay)
ksdelta=ksrenew(:,jj);
ksdelta(kk)=ksdelta(kk)+delta;
for i=1:Nfobs
[~,comres1(i,:),~]=fresolution(Nlay,2*pi*f(numbf(i)),exp(ksdelta),n_0,alpha,q_s,qp(numbf(i)),L);%ѡȡƵ�ʽ��м���
end
Amm1=angle(comres1);
virtAmap=zeros(Nfobs,Ndep);
for k=1:Ndep 
virtAmap(:,k)=Amm1(:,numdep(k));
end
virtapdelta=zeros(Nfobs,Ndep);
for i=1:(Nfobs)
for    j=1:Ndep
if(logth(numbf(i),numdep(j)))
 virtapdelta(i,j)=virtAmap(i,j);
end
end
end
virtapdelta=virtapdelta';
virtapdelta=virtapdelta(:);
virtapdelta(all(virtapdelta==0,2),:)=[];
jocob(:,kkk)=(virtapdelta-virta)/delta;
kkk=kkk+1;
end
%%%%%%%%%%%%%%
if(jj==1)
    for k=1:size(num_pam,2)
        for j=1:size(num_pam,2)
    Cyy(j,k)=0.66*0.66/2*exp(-(abs((num_pam(j)-num_pam(k)))*10/Lx).^2);
        end
    end
end
    
cdd=jocob*Cyy*jocob';
cdf=jocob*Cyy;
covdd=cdd+lambda(jj)*max(max(diag(cdd)))*eye(Ndd);
weigh=lsqminnorm(covdd,eye(Ndd))*cdf;
detal_par=weigh'*(refA-virta);


ksrenew(:,jj+1)=ksrenew(:,jj)+detal_par;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:Nfobs
[~,comres2(i,:),~]=fresolution(Nlay,2*pi*f(numbf(i)),exp(ksrenew(:,jj+1)),n_0,alpha,q_s,qp(numbf(i)),L);%ѡȡƵ�ʽ��м���
end
Am2=angle(comres2);
virtAma2=zeros(Nfobs,Ndep);
for k=1:Ndep 
virtAma2(:,k)=Am2(:,numdep(k));
end
virta2=zeros(Nfobs,Ndep);
for i=1:(Nfobs)
for    j=1:Ndep
if(logth(numbf(i),numdep(j)))
 virta2(i,j)=virtAma2(i,j);
end
end
end
virta2=virta2';
virta2=virta2(:);
virta2(all(virta2==0,2),:)=[];
rmse(jj+1)=sqrt(sum((virta2-refA).^2)/(Ndd));
if (rmse(jj+1)<rmse(jj))
lambda(jj+1)=lambda(jj)/lamb_adj;
Cyy=Cyy-weigh'*cdf;
else
    ksrenew(:,jj+1)=ksrenew(:,jj);
    lambda(jj+1)=lambda(jj)*(lamb_adj);
end
end

ksrmse=sqrt(sum((exp(ksrenew(:,Njtime+1))-K_s').^2)/19);
