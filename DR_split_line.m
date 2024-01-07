%% determine where DR=20, which means that on the bottom right of this line, the observation contains no information
Nlay=20;
K_s1=0.31*ones(1,Nlay); %Determine the K_s value for the soil 
alpha=11.35*ones(1,Nlay); %Determine the alpha value for the soil
n_0=0.39*ones(1,Nlay);  %Determine tehe porosity for the soil
L=[0.10,0.20,0.30,0.40,0.50,0.60,0.70,0.80,0.90,1.00,1.10,1.20,1.30,1.40,1.50,1.60,1.70,1.80,1.90,2.00];
deep=linspace(0.00,L(Nlay),L(Nlay)/0.01+1);
q_s=1/2*K_s1(1);
q_p=1/4*K_s1(1);
lgtt=linspace(-3,1,1001);
f=q_switch(K_s1(1),alpha(1),n_0(1),lgtt);
NN=length(f);
for j=1:NN
qp(j)=complex(q_p*cos(0),q_p*sin(0));
end 
%% 
N=1000;
t=linspace(0,N-1,N);
moni=zeros(N,1);
for i=1:NN
% [~,theta_p(i,:),~]=fresolution(Nlay,2*pi()*f(i),K_s,n_0,alpha,q_s,qp(i),L);
[~,theta_pA(i,:)]=fresolution(Nlay,2*pi()*f(i),K_s1,n_0,alpha,q_s,qp(i),L);

end


logth=log10(abs(theta_pA));
logth(find(logth(:)<(logth(1,1))-log10(20)))=0;
logth=logth';
% save('true.mat','Amm0','angm0','logth')
% subplot(2,1,1)
%  h=pcolor(res');
% set(h, 'LineStyle','none');
% colorbar;
% subplot(2,1,2)
% plot(moni,'LineWidth',3)
  currentFile =sprintf('logth.txt');
  save(currentFile,'logth','-ascii');
pp=lgtt;
  pp(2,:)=200;
  for i=1:1001
      for j=1:200
    if(logth(j,i)==0) 
  pp(2,i)=j-1;
  break;
    end
      end
  end
  pp=pp';
  currentFile =sprintf('pp.txt');
  save(currentFile,'pp','-ascii');