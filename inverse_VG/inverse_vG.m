%% load soil moisture data calculated by HYDRUS-1D equipped with VG model
thvg1=importdata('th1.csv');
thvg2=importdata('th2.csv');
thvg3=importdata('th15.csv');
thvg4=importdata('th30.csv');
thvg5=importdata('th90.csv');
thvg6=importdata('th180.csv');
thvg7=importdata('th365.csv');
T=[1,2,15,30,90,180,365]; %periods/d
tt=cell(1,7);
tt{1}=linspace(0,2400,2401)
tt{2}=linspace(0,2400,2401)
tt{3}=linspace(0,2400,2401)
tt{4}=linspace(0,2880,2881)
tt{5}=linspace(0,8640,8641)
tt{6}=linspace(0,12960,12961)
tt{7}=linspace(0,26280,26281)
for j=1:7
     eval(['thetavg=thvg',num2str(j),'.data']);
 time_array=1*T(j)*24+1:3*T(j)*24;
for i=1:201
zz=thetavg(i,time_array);
 t=tt{j}/24;
 [f,ff,fff]=fit(t(time_array)',(zz-mean(zz))','sin1');%sin function fitting
amp(j,i)=f.a1; %extract amplitude
 w(j,i)=f.b1;%extract period
ps(j,i)=f.c1;%extract phase shift
theta_s(j,i)=mean(zz);
end
end

Amm0=amp;
angm0=ps;
logth=ones(7,201);
save('true.mat','Amm0','angm0','theta_s','logth')


%% inverse VG
run colors_definitions
Nlay=20;
q_s=0.31/2;
q_p=0.31/4;
alpha=11.35*ones(1,Nlay); %alpha for GK model
n_0=0.39*ones(1,Nlay);
%Ks for VG model
K_s=[0.680161351572381;0.497247086918740;0.373301062344776;0.279876503832682;0.214912090075509;0.191280874090911;0.204443236385558;0.266741992554213;0.374640576891529;0.394988349720124;0.328978950255983;0.290465689688100;0.303832354535604;0.341712755916429;0.322657194749760;0.254185582915127;0.215177522155797;0.209692334891571;0.221149517873249;0.234554972631958];
Ndep=20;
numdep=linspace(5,195,20);
Njtime1=15;
T=[1,2,15,30,90,180,365]; %periods/d
f=1./T;
numbf=[1,4,6];
Nfobs=3;
ksrmse=zeros(4,4);
Lx=40;
K_s1=5*0.31;
ks=log(K_s1(1))*ones(Nlay,Njtime1+1);%Initial value
[ks(:,1:Njtime1+1),~]=inverse_amp(exp(ks(:,1)),Nfobs,numbf,Ndep,numdep,Njtime1,Lx,K_s,q_s,q_p,alpha,n_0);
Ksestima1(:,1)=exp(ks(:,Njtime1+1));
xlab=linspace(5,195,20);
L=[0.10,0.20,0.30,0.40,0.50,0.60,0.70,0.80,0.90,1.00,1.10,1.20,1.30,1.40,1.50,1.60,1.70,1.80,1.90,2.00];

fig=tiledlayout(2,4);
nexttile
for j=1:7
 for i=1:16
[~,vvv(i,:,j)]=fresolution(Nlay,2*pi*f(j),exp(ks(:,i)),n_0,alpha,q_s,complex(q_p*cos(0),q_p*sin(0)),L);
 end
end
for i=1:16
AA=reshape(vvv(i,numdep,numbf),[20,3]);
resm(i)=sqrt(sum((abs(AA')-Amm0(numbf,numdep)).^2,"all")/60);
end

plot(log10(resm),'*-','linewidth',2, 'color', color_scheme_aaas(10,:))
a=get(gca);
xmax=a.XLim;%
ymax=a.YLim;%
 set(gca,'FontSize',20,'FontWeight','bold');
    ylabel('\boldmath{$log(RMSE)$}','FontSize',20,'interpreter','latex')
   xlabel('\boldmath{iteration times $r$}','FontSize',20,'interpreter','latex')
 text(xmax(1)+0.01*(xmax(2)-xmax(1)),ymax(1)+0.95*(ymax(2)-ymax(1)),'(a)','FontSize',25,'FontWeight','bold');

nexttile
p1=plot(xlab,(K_s),'*-','linewidth',2, 'color', color_scheme_aaas(2,:));
hold on;
p2=plot(xlab,(Ksestima1),'*-','linewidth',2, 'color', color_scheme_aaas(10,:));
 set(gca,'FontSize',20,'FontWeight','bold');
 h1=legend([p1,p2],'{$K_{s}^{VG}$}',' {$K_{s}^{GK}$}','FontName','Times New Roman','interpreter','latex');
  set(h1,'Box','off')
  ylim([0 3.2])
a=get(gca);
xmax=a.XLim;%
ymax=a.YLim;%

text(xmax(1)-0.07*(xmax(2)-xmax(1)),ymax(1)+0*(ymax(2)-ymax(1)),'(b)','FontSize',25,'FontWeight','bold');
  ylabel('\boldmath{$K_{s}$}','FontSize',20,'interpreter','latex')
   xlabel('\boldmath{$z(cm)$}','FontSize',20,'interpreter','latex')
view(90,90)

nexttile([2 1])
p1=plot(linspace(0,200,201),abs(vvv(16,:,2)),'-','linewidth',2, 'color', color_scheme_aaas(8,:))
hold on;
color1=color_scheme_aaas(8,:);
sz=15
s1=scatter(linspace(0,200,201),Amm0(2,:),sz,color1,'filled');
p2=plot(linspace(0,200,201),abs(vvv(16,:,7)),'-','linewidth',2, 'color', color_scheme_aaas(7,:))
hold on;
color2=color_scheme_aaas(7,:);
s2=scatter(linspace(0,200,201),Amm0(7,:),sz,color2,'filled');
 set(gca,'FontSize',20,'FontWeight','bold');

  ylim([0.0128 0.016])
a=get(gca);
xmax=a.XLim;%
ymax=a.YLim;%
text(xmax(1)-0.025*(xmax(2)-xmax(1)),ymax(1)+0*(ymax(2)-ymax(1)),'(d)','FontSize',25,'FontWeight','bold');
   ylabel('\boldmath{$Amp(\theta$)}','FontSize',20,'interpreter','latex')
   xlabel('\boldmath{$z(cm)$}','FontSize',20,'interpreter','latex')
view(90,90)

nexttile([2 1])
p1=plot(linspace(0,200,201),angle(vvv(16,:,2)),'-','linewidth',2, 'color', color_scheme_aaas(8,:))
hold on;
color1=color_scheme_aaas(8,:);
sz=15
s1=scatter(linspace(0,200,201),angm0(2,:),sz,color1,'filled');
p2=plot(linspace(0,200,201),angle(vvv(16,:,7)),'-','linewidth',2, 'color', color_scheme_aaas(7,:))
hold on;
color2=color_scheme_aaas(7,:);
s2=scatter(linspace(0,200,201),angm0(7,:),sz,color2,'filled');
 set(gca,'FontSize',20,'FontWeight','bold');
 h1=legend([p1,s1,p2,s2],'{Simulated by GK model for $t_{p}=2d$}',' {Extracted by VG model for $t_{p}=2d$}','{Simulated by GK model for $t_{p}=365d$}',' {Extracted by VG model for $t_{p}=365d$}','FontName','Times New Roman','interpreter','latex');
  set(h1,'Box','off')
  set(h1,'FontSize',15)
h1.Location='north';
%   ylim([0.0128 0.018])
a=get(gca);
xmax=a.XLim;%获取横坐标上下限
ymax=a.YLim;%获取纵坐标上下限
text(xmax(1)-0.025*(xmax(2)-xmax(1)),ymax(1)+0*(ymax(2)-ymax(1)),'(e)','FontSize',25,'FontWeight','bold');
 ylabel('\boldmath{$PS(\theta$)}','FontSize',20,'interpreter','latex')
   xlabel('\boldmath{$z(cm)$}','FontSize',20,'interpreter','latex')
view(90,90)


nexttile([1 2])
x=Ksestima1;
y=K_s;
mdl = fitlm(K_s,Ksestima1);
p=polyfit(x,y,1);%
yfit=polyval(p,x);%
a = num2str(p(1),'%.3f');%
b = num2str(p(2),'%.3f');%
r2 = num2str(mdl.Rsquared.Ordinary,'%.3f')
% mdl.Coefficients
anova(mdl,'summary')
pp=plot(mdl);
set(gca,'FontSize',20,'FontWeight','bold');
Formu = ['y=',a,'x',b,'  R^2=',r2];%
xlim([0.15 0.7])
ylim([1.9 2.5])
text(0.5,2.1,Formu,'FontSize',20,'FontWeight','bold','color',color_scheme_npg(4,:)) 
pp(1).LineWidth=3;
pp(2).LineWidth=3;
pp(3).LineWidth=3;
pp(4).LineWidth=3;
pp(1).Color=color_scheme_npg(3,:);
pp(2).Color=color_scheme_npg(4,:);
pp(3).Color=color_scheme_npg(6,:);
pp(4).Color=color_scheme_npg(6,:);
  ylabel('\boldmath{$K_{s}^{GK}$}','FontSize',20,'interpreter','latex')
   xlabel('\boldmath{$K_{s}^{VG}$}','FontSize',20,'interpreter','latex')
h=legend('data','Fit','95% conf.bounds');
  set(h,'Box','off')
  set(h,'FontSize',20)
  h.Location='northwest';
  a=get(gca);
xmax=a.XLim;%
ymax=a.YLim;%
 text(xmax(1)+0.01*(xmax(2)-xmax(1)),ymax(1)+1.05*(ymax(2)-ymax(1)),'(c)','FontSize',25,'FontWeight','bold');
title(' ','FontSize',20,'FontWeight','bold')












