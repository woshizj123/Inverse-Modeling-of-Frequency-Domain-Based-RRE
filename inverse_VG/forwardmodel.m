%% ����ģ��
% ��������ģ�ͣ�ȷ��ʵ��ֵ��ͬ���ʵ��ֵ�Ĵ�С�������ֲ�0-10cm��10-30cm,30-85cm,85-150cm�����������۲�ֲ�0-10cm������20cmһ���ֲ㣬ģ��150cm
% �����۲�ʵ����5cm,15cm,45cm,95cm����������ÿ1cm���ù۲�
%% ���䴢����,��������
function [ ]=forwardmodel(K_s,K_s1,alpha,n_0,q_s,q_p)
Nlay=20;%��������
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
%% ����ģ��
N=1000;
t=linspace(0,N-1,N);
moni=zeros(N,1);
for i=1:NN
[~,theta_p(i,:),~]=fresolution(Nlay,2*pi()*f(i),K_s,n_0,alpha,q_s,qp(i),L);%����ֱ��Ӧ˥����ʵ������ÿ����̬��Ⱥ�ˮ�ʣ�������ˮ�ʣ�����������λ�����ĳһ���ʱ�䣻����˥�����ӣ�����ȵݼ���
[~,theta_pA(i,:),~]=fresolution(Nlay,2*pi()*f(i),K_s1,n_0,alpha,q_s,qp(i),L);%����ֱ��Ӧ˥����ʵ������ÿ����̬��Ⱥ�ˮ�ʣ�������ˮ�ʣ�����������λ�����ĳһ���ʱ�䣻����˥�����ӣ�����ȵݼ���

end
[theta_s]=[]%����ֱ��Ӧ˥����ʵ������ÿ����̬��Ⱥ�ˮ�ʣ�������ˮ�ʣ�����������λ�����ĳһ���ʱ�䣻����˥�����ӣ�����ȵݼ���
res=repmat(theta_s,N,1);

% for i=1:NN
%     moni(:,1)=moni(:,1)+abs(theta_p0(i))*cos(2*pi*f(i)*t'+angle(theta_p0(i)));
% end

Amm0=abs(theta_p);%����Ƶ�ʶ�Ӧ�����
angm0=angle(theta_p);%����Ƶ���¶�Ӧ�ķ���
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