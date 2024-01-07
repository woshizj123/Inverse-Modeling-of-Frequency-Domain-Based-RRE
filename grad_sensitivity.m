%% parameter determines
Nlay=200;%soil layers
K_s=0.31*ones(1,Nlay); %Determine the Ks value for every layers
alpha=11.35*ones(1,Nlay); %Determine the alpha value for every layers
 h_e=-0.244*ones(1,Nlay);%Determine the he for every layers
n_0=0.39*ones(1,Nlay);  %Determine tehe porosity for every layers
% mu=2.719*ones(1,Nlay); %Determine the fitting parameter for every layers 
L=linspace(0.01,2.00,Nlay); %grids
% deep=linspace(0.1,L(Nlay),L(Nlay)/0.1);
q_s=1/2*K_s(1); %steady flux
q_p=1/4*K_s(1); %fluctuation flux
lgtt=linspace(-3,1,1001); %Townley number
f=q_switch(K_s(1),alpha(1),n_0(1),lgtt); % Townley number Convert to frequency f
Nf=size(f,2);
%% 
%%%%% jocob matrix
for j=1:Nf
for i=1:200
    K_s=0.31*ones(1,Nlay); 
    [~,dampthetap]=fresolution(Nlay,2*pi()*f(j),K_s,n_0,alpha,q_s,q_p,L);%
    K_s(i)=K_s(i)+0.0001;
[~,dampthetap1]=fresolution(Nlay,2*pi()*f(j),K_s,n_0,alpha,q_s,q_p,L);%
aa1(i,j)=(abs(dampthetap1(1))-abs(dampthetap(1)))/0.0001;
bb1(i,j)=(angle(dampthetap1(1))-angle(dampthetap(1)))/0.0001;
aa2(i,j)=(abs(dampthetap1(100))-abs(dampthetap(100)))/0.0001;
bb2(i,j)=(angle(dampthetap1(100))-angle(dampthetap(100)))/0.0001;
aa3(i,j)=(abs(dampthetap1(180))-abs(dampthetap(180)))/0.0001;
bb3(i,j)=(angle(dampthetap1(180))-angle(dampthetap(180)))/0.0001;
end
end
for i=1:1001
    bb1(:,i)=bb1(:,i)/f(i);
    bb2(:,i)=bb2(:,i)/f(i);
    bb3(:,i)=bb3(:,i)/f(i);
end
%%%%%%%%%% Cdy


Lx=40;
    for k=1:Nlay
        for j=1:Nlay
    Cyy(j,k)=0.66*0.66*exp(-(abs(((j)-(k)))*10/Lx).^2); % Parameter variance
        end
    end
    
    for k=1:3
   eval(['bot',num2str(k),'=','transpose(aa',num2str(k),')','*Cyy*aa',num2str(k)]);
   eval(['top',num2str(k),'=','Cyy*aa',num2str(k)]);
   eval(['angbot',num2str(k),'=','transpose(bb',num2str(k),')','*Cyy*bb',num2str(k)]);
   eval(['angtop',num2str(k),'=','Cyy*bb',num2str(k)]);
    end
    
 %% cross-correlation matrix ρ
  for i=1:200
        for j=1:1001
          amp1(i,j)=top1(i,j)/sqrt(bot1(j,j)*Cyy(i,i));
                    amp2(i,j)=top2(i,j)/sqrt(bot2(j,j)*Cyy(i,i));
                              amp3(i,j)=top3(i,j)/sqrt(bot3(j,j)*Cyy(i,i));
           ang1(i,j)=angtop1(i,j)/sqrt(angbot1(j,j)*Cyy(i,i));
                    ang2(i,j)=angtop2(i,j)/sqrt(angbot2(j,j)*Cyy(i,i));
                              ang3(i,j)=angtop3(i,j)/sqrt(angbot3(j,j)*Cyy(i,i));

        end
  end
       
%   save('cor_Am.mat','amp1','amp2','amp3')
% save('cor_Ang.mat','ang1','ang2','ang3')
%  currentFile =sprintf('amp1.txt');
%  save(currentFile,'amp1','-ascii');
%  currentFile =sprintf('amp2.txt');
%  save(currentFile,'amp2','-ascii');
%  currentFile =sprintf('amp3.txt');
%  save(currentFile,'amp3','-ascii');
%   currentFile =sprintf('ang1.txt');
%  save(currentFile,'ang1','-ascii');
%  currentFile =sprintf('ang2.txt');
%  save(currentFile,'ang2','-ascii');
%  currentFile =sprintf('ang3.txt');
%  save(currentFile,'ang3','-ascii');          
            
            
