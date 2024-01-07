%% 正演模型
% 利用正演模型，确定实际值不同深度实际值的大小，土壤分层0-10cm，10-30cm,30-85cm,85-150cm，虚拟算例观测分层0-10cm，后面20cm一个分层，模拟150cm
% 土壤观测实际在5cm,15cm,45cm,95cm处，虚拟中每1cm设置观测
%% 分配储存区,参数设置
function [ ]=forwardmodel(K_s,K_s1,alpha,n_0,q_s,q_p)
Nlay=20;%土壤层数
L=[0.10,0.20,0.30,0.40,0.50,0.60,0.70,0.80,0.90,1.00,1.10,1.20,1.30,1.40,1.50,1.60,1.70,1.80,1.90,2.00];
deep=linspace(0.0,L(Nlay),L(Nlay)/0.01+1);
% q_s=5.e-03;
% q_p=1.e-04;
lgtt=linspace(-3,1,1001);
f=q_switch(K_s1(1),alpha(1),n_0(1),lgtt);
NN=length(f);
for j=1:NN
qp(j)=complex(q_p*cos(0),q_p*sin(0));
end 
%% 正演模型
N=1000;
t=linspace(0,N-1,N);
moni=zeros(N,1);
for i=1:NN
[~,theta_p(i,:),~]=fresolution(Nlay,2*pi()*f(i),K_s,n_0,alpha,q_s,qp(i),L);%输出分别对应衰减（实数）；每个稳态深度含水率；波动含水率（复数）；相位差；到达某一深度时间；波动衰减因子（随深度递减）
[~,theta_pA(i,:),~]=fresolution(Nlay,2*pi()*f(i),K_s1,n_0,alpha,q_s,qp(i),L);%输出分别对应衰减（实数）；每个稳态深度含水率；波动含水率（复数）；相位差；到达某一深度时间；波动衰减因子（随深度递减）

end
[theta_s]=[]%输出分别对应衰减（实数）；每个稳态深度含水率；波动含水率（复数）；相位差；到达某一深度时间；波动衰减因子（随深度递减）
res=repmat(theta_s,N,1);

% for i=1:NN
%     moni(:,1)=moni(:,1)+abs(theta_p0(i))*cos(2*pi*f(i)*t'+angle(theta_p0(i)));
% end

Amm0=abs(theta_p);%所有频率对应的振幅
angm0=angle(theta_p);%所有频率下对应的幅角
% for k=1:L(Nlay)/0.01
% for i=1:NN
%      for j=2:Nlay
%     if deep(k)<=L(j)&&deep(k)>L(j-1)
%     res(:,k)=res(:,k)+Amm0(i,k)*cos(2*pi*f(i)*t'+angm0(i,k));
%     end
%     end
% end
% end
% for k=1:L(1)/0.01
%     for i=1:NN
% res(:,k)=res(:,k)+Amm0(i,k)*cos(2*pi*f(i)*t'+angm0(i,k));
%     end
% end

logth=log10(abs(theta_pA));
logth(find(logth(:)<(logth(1,1))-1.3))=0;
save('true.mat','Amm0','angm0','theta_s','logth')
% subplot(2,1,1)
%  h=pcolor(res');
% set(h, 'LineStyle','none');
% colorbar;
% subplot(2,1,2)
% plot(moni,'LineWidth',3)
%  currentFile =sprintf('logth.txt');
%  save(currentFile,'logth','-ascii');