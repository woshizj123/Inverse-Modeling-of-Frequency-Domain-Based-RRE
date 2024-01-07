%% Effects of different observation scheme settings (different numbers of spatial observations or frequenices)
texture=readmatrix('texture.xlsx','Sheet','Sheet1','Range','B2:G6');
Nlay=20;
alpha=texture(4,1)*ones(1,Nlay); %Determine the alpha value for the soil
n_0=texture(4,2)*ones(1,Nlay);  %Determine tehe porosity for the soil
L=[0.10,0.20,0.30,0.40,0.50,0.60,0.70,0.80,0.90,1.00,1.10,1.20,1.30,1.40,1.50,1.60,1.70,1.80,1.90,2.00];
K_s1=texture(4,4)*ones(1,Nlay); %Determine the K_s value for the soil 
kk=[10,20,30,40,50,80,100];
q_s=1/2*K_s1(1);
q_p=1/4*K_s1(1);
[Nfobs1,numbf1]=fre_arry2(41);
lgtt=linspace(-3,1,1001);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%nd=10,nf=1
tt=tiledlayout(3,3);
nexttile
K_s=[0.680161351572381;0.497247086918740;0.373301062344776;0.279876503832682;0.214912090075509;0.191280874090911;0.204443236385558;0.266741992554213;0.374640576891529;0.394988349720124;0.328978950255983;0.290465689688100;0.303832354535604;0.341712755916429;0.322657194749760;0.254185582915127;0.215177522155797;0.209692334891571;0.221149517873249;0.234554972631958];
for j=1:1
    numbf=numbf1(1:1);
    Nfobs=length(numbf);
  a_lookf=lgtt(numbf);
%Ndep=5;
Ndep=10;
numdep=linspace(5,185,10);
Njtime1=10;
Njtime2=20;
Njtime3=20;
ksrmse=zeros(4,4);
Lx=40;
ks=K_s1(1)*ones(Nlay,Njtime2+Njtime3+Njtime1+3);
forwardmodel(K_s(:,j),K_s1,alpha,n_0,q_s,q_p);
[ks(:,1:Njtime1+1),~]=inverse_steady(ks(:,1),1,numbf,Ndep,numdep,Njtime1,Lx,K_s(:,j),q_s,q_p,alpha,n_0);
[ks(:,Njtime1+2:Njtime2+Njtime1+2),~]=inverse_amp(exp(ks(:,Njtime1+1)),Nfobs,numbf,Ndep,numdep,Njtime2,Lx,K_s(:,j),q_s,q_p,alpha,n_0);
[ks(:,Njtime2+Njtime1+3:Njtime3+Njtime2+Njtime1+3),~]=inverse_ps(exp(ks(:,Njtime2+Njtime1+2)),Nfobs,numbf,Ndep,numdep,Njtime3,Lx,K_s(:,j),q_s,q_p,alpha,n_0);
Ksestima1(:,j)=exp(ks(:,Njtime1+1));
Ksestima2(:,j)=exp(ks(:,Njtime2+Njtime1+2));
Ksestima3(:,j)=exp(ks(:,Njtime2+Njtime3+Njtime1+3));
end
run colors_definitions
rmse0=abs(sqrt(sum((log(0.31)-log(K_s)).^2)/19)./mean(log(K_s)));
rmse1=abs(sqrt(sum((log(Ksestima1)-log(K_s)).^2)/19)./mean(log(K_s)));
rmse2=abs(sqrt(sum((log(Ksestima2)-log(K_s)).^2)/19)./mean(log(K_s)));
rmse3=abs(sqrt(sum((log(Ksestima3)-log(K_s)).^2)/19)./mean(log(K_s)));
xlab=linspace(5,195,20);
p1=plot(xlab,log(K_s),'*-','linewidth',2, 'color', color_scheme_aaas(10,:));
hold on;
p2=plot(xlab,log(Ksestima1),'*-','linewidth',2, 'color', color_scheme_aaas(2,:));
p3=plot(xlab,log(Ksestima2),'*-','linewidth',2, 'color', color_scheme_aaas(3,:));
p4=plot(xlab,log(Ksestima3),'*-','linewidth',2, 'color', color_scheme_aaas(4,:));
ylim([-1.8 -0.2])
  a=get(gca);
  xmax=a.XLim;
  ymax=a.YLim;
  p5=scatter(numdep,ymax(1),'m');

  set(gca,'FontSize',15,'FontWeight','bold');
   set(gca,'box','off');

   str1={['$$NRMSE_1=$$' num2str(rmse1,'%4.4f') ]};
  str2={['$$NRMSE_2=$$' num2str(rmse2,'%4.4f') ]};
   str3={['$$NRMSE_3=$$' num2str(rmse3,'%4.4f') ]};
   str4={['$$n_f=$$',num2str(Nfobs),',$$n_{sp}=$$',num2str(Ndep)]};

   kk1=[0.15 0.7];
  x0=xmax(1)+kk1(1)*(xmax(2)-xmax(1));
  y0=ymax(1)+kk1(2)*(ymax(2)-ymax(1));
   kk2=[0.25 0.7];
  x1=xmax(1)+kk2(1)*(xmax(2)-xmax(1));
  y1=ymax(1)+kk2(2)*(ymax(2)-ymax(1));
     kk3=[0.35 0.7];
  x2=xmax(1)+kk3(1)*(xmax(2)-xmax(1));
  y2=ymax(1)+kk3(2)*(ymax(2)-ymax(1));
  text(x0,y0,str1,'FontSize',12,'FontWeight','bold', 'color', color_scheme_aaas(2,:),'interpreter','latex');
  text(x1,y1,str2,'FontSize',12,'FontWeight','bold', 'color', color_scheme_aaas(3,:),'interpreter','latex');
  text(x2,y2,str3,'FontSize',12,'FontWeight','bold', 'color', color_scheme_aaas(4,:),'interpreter','latex');
view(90,90)
subtitle(str4,'FontSize',25,'FontWeight','bold','interpreter','latex')
   title('(a)','position',[xmax(1)-8,ymax(1)],'FontSize',20);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%nd=10,nf=1

nexttile
for j=1:1
    numbf=numbf1(25:25);
    Nfobs=length(numbf);
  a_lookf=lgtt(numbf);
Ndep=10;
numdep=linspace(5,185,10);
Njtime1=10;
Njtime2=20;
Njtime3=20;
ksrmse=zeros(4,4);
Lx=40;
ks=K_s1(1)*ones(Nlay,Njtime2+Njtime3+Njtime1+3);
forwardmodel(K_s(:,j),K_s1,alpha,n_0,q_s,q_p);
[ks(:,1:Njtime1+1),~]=inverse_steady(ks(:,1),1,numbf,Ndep,numdep,Njtime1,Lx,K_s(:,j),q_s,q_p,alpha,n_0);
[ks(:,Njtime1+2:Njtime2+Njtime1+2),~]=inverse_amp(exp(ks(:,Njtime1+1)),Nfobs,numbf,Ndep,numdep,Njtime2,Lx,K_s(:,j),q_s,q_p,alpha,n_0);
[ks(:,Njtime2+Njtime1+3:Njtime3+Njtime2+Njtime1+3),~]=inverse_ps(exp(ks(:,Njtime2+Njtime1+2)),Nfobs,numbf,Ndep,numdep,Njtime3,Lx,K_s(:,j),q_s,q_p,alpha,n_0);
Ksestima1(:,j)=exp(ks(:,Njtime1+1));
Ksestima2(:,j)=exp(ks(:,Njtime2+Njtime1+2));
Ksestima3(:,j)=exp(ks(:,Njtime2+Njtime3+Njtime1+3));
end
rmse0=abs(sqrt(sum((log(0.31)-log(K_s)).^2)/19)./mean(log(K_s)));
rmse1=abs(sqrt(sum((log(Ksestima1)-log(K_s)).^2)/19)./mean(log(K_s)));
rmse2=abs(sqrt(sum((log(Ksestima2)-log(K_s)).^2)/19)./mean(log(K_s)));
rmse3=abs(sqrt(sum((log(Ksestima3)-log(K_s)).^2)/19)./mean(log(K_s)));

xlab=linspace(5,195,20);
p1=plot(xlab,log(K_s),'*-','linewidth',2, 'color', color_scheme_aaas(10,:));
hold on;
p2=plot(xlab,log(Ksestima1),'*-','linewidth',2, 'color', color_scheme_aaas(2,:));
p3=plot(xlab,log(Ksestima2),'*-','linewidth',2, 'color', color_scheme_aaas(3,:));
p4=plot(xlab,log(Ksestima3),'*-','linewidth',2, 'color', color_scheme_aaas(4,:));
ylim([-1.8 -0.2])
  a=get(gca);
  xmax=a.XLim;
  ymax=a.YLim;
  set(gca,'FontSize',15,'FontWeight','bold');
p5=scatter(numdep,ymax(1),'m');
    set(gca,'box','off');
    set(gca,'xticklabel',[])

   str1={['$$NRMSE_1=$$' num2str(rmse1,'%4.4f') ]};
  str2={['$$NRMSE_2=$$' num2str(rmse2,'%4.4f') ]};
   str3={['$$NRMSE_3=$$' num2str(rmse3,'%4.4f') ]};
   str4={['$$n_f=$$',num2str(Nfobs),',$$n_{sp}=$$',num2str(Ndep)]};
   kk1=[0.15 0.7];
  x0=xmax(1)+kk1(1)*(xmax(2)-xmax(1));
  y0=ymax(1)+kk1(2)*(ymax(2)-ymax(1));
   kk2=[0.25 0.7];
  x1=xmax(1)+kk2(1)*(xmax(2)-xmax(1));
  y1=ymax(1)+kk2(2)*(ymax(2)-ymax(1));
     kk3=[0.35 0.7];
  x2=xmax(1)+kk3(1)*(xmax(2)-xmax(1));
  y2=ymax(1)+kk3(2)*(ymax(2)-ymax(1));
  text(x0,y0,str1,'FontSize',12,'FontWeight','bold', 'color', color_scheme_aaas(2,:),'interpreter','latex');
  text(x1,y1,str2,'FontSize',12,'FontWeight','bold', 'color', color_scheme_aaas(3,:),'interpreter','latex');
  text(x2,y2,str3,'FontSize',12,'FontWeight','bold', 'color', color_scheme_aaas(4,:),'interpreter','latex');
view(90,90)
subtitle(str4,'FontSize',25,'FontWeight','bold','interpreter','latex')
   title('(b)','position',[xmax(1)-8,ymax(1)],'FontSize',20);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%nd=10,nf=2
nexttile
for j=1:1
    numbf=numbf1([1 25]);
    Nfobs=length(numbf);
  a_lookf=lgtt(numbf);
Ndep=10;
numdep=linspace(5,185,10);
%numdep=linspace(5,195,20);
%numdep=[5];
Njtime1=10;
Njtime2=20;
Njtime3=20;
ksrmse=zeros(4,4);
Lx=40;
ks=K_s1(1)*ones(Nlay,Njtime2+Njtime3+Njtime1+3);
forwardmodel(K_s(:,j),K_s1,alpha,n_0,q_s,q_p);
[ks(:,1:Njtime1+1),~]=inverse_steady(ks(:,1),1,numbf,Ndep,numdep,Njtime1,Lx,K_s(:,j),q_s,q_p,alpha,n_0);
[ks(:,Njtime1+2:Njtime2+Njtime1+2),~]=inverse_amp(exp(ks(:,Njtime1+1)),Nfobs,numbf,Ndep,numdep,Njtime2,Lx,K_s(:,j),q_s,q_p,alpha,n_0);
[ks(:,Njtime2+Njtime1+3:Njtime3+Njtime2+Njtime1+3),~]=inverse_ps(exp(ks(:,Njtime2+Njtime1+2)),Nfobs,numbf,Ndep,numdep,Njtime3,Lx,K_s(:,j),q_s,q_p,alpha,n_0);
Ksestima1(:,j)=exp(ks(:,Njtime1+1));
Ksestima2(:,j)=exp(ks(:,Njtime2+Njtime1+2));
Ksestima3(:,j)=exp(ks(:,Njtime2+Njtime3+Njtime1+3));
end
rmse0=abs(sqrt(sum((log(0.31)-log(K_s)).^2)/19)./mean(log(K_s)));
rmse1=abs(sqrt(sum((log(Ksestima1)-log(K_s)).^2)/19)./mean(log(K_s)));
rmse2=abs(sqrt(sum((log(Ksestima2)-log(K_s)).^2)/19)./mean(log(K_s)));
rmse3=abs(sqrt(sum((log(Ksestima3)-log(K_s)).^2)/19)./mean(log(K_s)));

xlab=linspace(5,195,20);
p1=plot(xlab,log(K_s),'*-','linewidth',2, 'color', color_scheme_aaas(10,:));
hold on;
p2=plot(xlab,log(Ksestima1),'*-','linewidth',2, 'color', color_scheme_aaas(2,:));
p3=plot(xlab,log(Ksestima2),'*-','linewidth',2, 'color', color_scheme_aaas(3,:));
p4=plot(xlab,log(Ksestima3),'*-','linewidth',2, 'color', color_scheme_aaas(4,:));
% ylim([-1.8 -0.2])
ylim([-1.8 -0.2])
   a=get(gca);
  xmax=a.XLim;
  ymax=a.YLim;
  set(gca,'FontSize',15,'FontWeight','bold');
p5=scatter(numdep,ymax(1),'m');
    set(gca,'box','off');
set(gca,'xticklabel',[])

   str1={['$$NRMSE_1=$$' num2str(rmse1,'%4.4f') ]};
  str2={['$$NRMSE_2=$$' num2str(rmse2,'%4.4f') ]};
   str3={['$$NRMSE_3=$$' num2str(rmse3,'%4.4f') ]};
   str4={['$$n_f=$$',num2str(Nfobs),',$$n_{sp}=$$',num2str(Ndep)]};

   kk1=[0.15 0.7];%
  x0=xmax(1)+kk1(1)*(xmax(2)-xmax(1));%
  y0=ymax(1)+kk1(2)*(ymax(2)-ymax(1));%
   kk2=[0.25 0.7];%
  x1=xmax(1)+kk2(1)*(xmax(2)-xmax(1));%
  y1=ymax(1)+kk2(2)*(ymax(2)-ymax(1));%
     kk3=[0.35 0.7];%
  x2=xmax(1)+kk3(1)*(xmax(2)-xmax(1));%
  y2=ymax(1)+kk3(2)*(ymax(2)-ymax(1));%
  text(x0,y0,str1,'FontSize',12,'FontWeight','bold', 'color', color_scheme_aaas(2,:),'interpreter','latex');
  text(x1,y1,str2,'FontSize',12,'FontWeight','bold', 'color', color_scheme_aaas(3,:),'interpreter','latex');
  text(x2,y2,str3,'FontSize',12,'FontWeight','bold', 'color', color_scheme_aaas(4,:),'interpreter','latex');
view(90,90)
subtitle(str4,'FontSize',25,'FontWeight','bold','interpreter','latex')
   title('(c)','position',[xmax(1)-8,ymax(1)],'FontSize',20);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%nd=1,nf=21
nexttile
for j=1:1
        numbf=linspace(1,501,21);
    Nfobs=length(numbf);
  a_lookf=lgtt(numbf);
%Ndep=5;
Ndep=1;
%numdep=linspace(5,195,20);
numdep=[1];
Njtime1=10;
Njtime2=20;
Njtime3=20;
ksrmse=zeros(4,4);
Lx=40;
ks=K_s1(1)*ones(Nlay,Njtime2+Njtime3+Njtime1+3);
forwardmodel(K_s(:,j),K_s1,alpha,n_0,q_s,q_p);
[ks(:,1:Njtime1+1),~]=inverse_steady(ks(:,1),1,numbf,Ndep,numdep,Njtime1,Lx,K_s(:,j),q_s,q_p,alpha,n_0);
[ks(:,Njtime1+2:Njtime2+Njtime1+2),~]=inverse_amp(exp(ks(:,Njtime1+1)),Nfobs,numbf,Ndep,numdep,Njtime2,Lx,K_s(:,j),q_s,q_p,alpha,n_0);
[ks(:,Njtime2+Njtime1+3:Njtime3+Njtime2+Njtime1+3),~]=inverse_ps(exp(ks(:,Njtime2+Njtime1+2)),Nfobs,numbf,Ndep,numdep,Njtime3,Lx,K_s(:,j),q_s,q_p,alpha,n_0);
Ksestima1(:,j)=exp(ks(:,Njtime1+1));
Ksestima2(:,j)=exp(ks(:,Njtime2+Njtime1+2));
Ksestima3(:,j)=exp(ks(:,Njtime2+Njtime3+Njtime1+3));
end
cmap = get(0, 'defaultaxescolororder');
rmse0=abs(sqrt(sum((log(0.31)-log(K_s)).^2)/19)./mean(log(K_s)));
rmse1=abs(sqrt(sum((log(Ksestima1)-log(K_s)).^2)/19)./mean(log(K_s)));
rmse2=abs(sqrt(sum((log(Ksestima2)-log(K_s)).^2)/19)./mean(log(K_s)));
rmse3=abs(sqrt(sum((log(Ksestima3)-log(K_s)).^2)/19)./mean(log(K_s)));
xlabel('$$z$$(cm)','FontSize',25,'FontWeight','bold','interpreter','latex');

xlab=linspace(5,195,20);
p1=plot(xlab,log(K_s),'*-','linewidth',2, 'color', color_scheme_aaas(10,:));
hold on;
p2=plot(xlab,log(Ksestima1),'*-','linewidth',2, 'color', color_scheme_aaas(2,:));
p3=plot(xlab,log(Ksestima2),'*-','linewidth',2, 'color', color_scheme_aaas(3,:));
p4=plot(xlab,log(Ksestima3),'*-','linewidth',2, 'color', color_scheme_aaas(4,:));
ylim([-2 0])
  a=get(gca);
  xmax=a.XLim;%
  ymax=a.YLim;%
  set(gca,'FontSize',15,'FontWeight','bold');
p5=scatter(numdep,ymax(1),'m');
    set(gca,'box','off');

   str1={['$$NRMSE_1=$$' num2str(rmse1,'%4.4f') ]};
  str2={['$$NRMSE_2=$$' num2str(rmse2,'%4.4f') ]};
   str3={['$$NRMSE_3=$$' num2str(rmse3,'%4.4f') ]};
   str4={['$$n_f=$$',num2str(Nfobs),',$$n_{sp}=$$',num2str(Ndep)]};

   kk1=[0.15 0.7];%
  x0=xmax(1)+kk1(1)*(xmax(2)-xmax(1));%
  y0=ymax(1)+kk1(2)*(ymax(2)-ymax(1));
   kk2=[0.25 0.7];%
  x1=xmax(1)+kk2(1)*(xmax(2)-xmax(1));%
  y1=ymax(1)+kk2(2)*(ymax(2)-ymax(1));%
     kk3=[0.35 0.7];%
  x2=xmax(1)+kk3(1)*(xmax(2)-xmax(1));%
  y2=ymax(1)+kk3(2)*(ymax(2)-ymax(1));%
  text(x0,y0,str1,'FontSize',12,'FontWeight','bold', 'color', color_scheme_aaas(2,:),'interpreter','latex');
  text(x1,y1,str2,'FontSize',12,'FontWeight','bold', 'color', color_scheme_aaas(3,:),'interpreter','latex');
  text(x2,y2,str3,'FontSize',12,'FontWeight','bold', 'color', color_scheme_aaas(4,:),'interpreter','latex');
view(90,90)
subtitle(str4,'FontSize',25,'FontWeight','bold','interpreter','latex')
   title('(d)','position',[xmax(1)-8,ymax(1)],'FontSize',20);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%nd=1,nf=21
nexttile
for j=1:1
        numbf=linspace(501,1001,21);
    Nfobs=length(numbf);
  a_lookf=lgtt(numbf);
%Ndep=5;
Ndep=1;
%numdep=linspace(5,195,20);
numdep=[1];
Njtime1=10;
Njtime2=20;
Njtime3=20;
ksrmse=zeros(4,4);
Lx=40;
ks=K_s1(1)*ones(Nlay,Njtime2+Njtime3+Njtime1+3);
forwardmodel(K_s(:,j),K_s1,alpha,n_0,q_s,q_p);
[ks(:,1:Njtime1+1),~]=inverse_steady(ks(:,1),1,numbf,Ndep,numdep,Njtime1,Lx,K_s(:,j),q_s,q_p,alpha,n_0);
[ks(:,Njtime1+2:Njtime2+Njtime1+2),~]=inverse_amp(exp(ks(:,Njtime1+1)),Nfobs,numbf,Ndep,numdep,Njtime2,Lx,K_s(:,j),q_s,q_p,alpha,n_0);
[ks(:,Njtime2+Njtime1+3:Njtime3+Njtime2+Njtime1+3),~]=inverse_ps(exp(ks(:,Njtime2+Njtime1+2)),Nfobs,numbf,Ndep,numdep,Njtime3,Lx,K_s(:,j),q_s,q_p,alpha,n_0);
Ksestima1(:,j)=exp(ks(:,Njtime1+1));
Ksestima2(:,j)=exp(ks(:,Njtime2+Njtime1+2));
Ksestima3(:,j)=exp(ks(:,Njtime2+Njtime3+Njtime1+3));
end
cmap = get(0, 'defaultaxescolororder');
rmse0=abs(sqrt(sum((log(0.31)-log(K_s)).^2)/19)./mean(log(K_s)));
rmse1=abs(sqrt(sum((log(Ksestima1)-log(K_s)).^2)/19)./mean(log(K_s)));
rmse2=abs(sqrt(sum((log(Ksestima2)-log(K_s)).^2)/19)./mean(log(K_s)));
rmse3=abs(sqrt(sum((log(Ksestima3)-log(K_s)).^2)/19)./mean(log(K_s)));

xlab=linspace(5,195,20);
p1=plot(xlab,log(K_s),'*-','linewidth',2, 'color', color_scheme_aaas(10,:));
hold on;
p2=plot(xlab,log(Ksestima1),'*-','linewidth',2, 'color', color_scheme_aaas(2,:));
p3=plot(xlab,log(Ksestima2),'*-','linewidth',2, 'color', color_scheme_aaas(3,:));
p4=plot(xlab,log(Ksestima3),'*-','linewidth',2, 'color', color_scheme_aaas(4,:));
ylim([-2 0])
  a=get(gca);
  xmax=a.XLim;%
  ymax=a.YLim;%
  set(gca,'FontSize',15,'FontWeight','bold');
p5=scatter(numdep,ymax(1),'m');
    set(gca,'box','off');
set(gca,'xticklabel',[])

   str1={['$$NRMSE_1=$$' num2str(rmse1,'%4.4f') ]};
  str2={['$$NRMSE_2=$$' num2str(rmse2,'%4.4f') ]};
   str3={['$$NRMSE_3=$$' num2str(rmse3,'%4.4f') ]};
   str4={['$$n_f=$$',num2str(Nfobs),',$$n_{sp}=$$',num2str(Ndep)]};

   kk1=[0.15 0.7];%
  x0=xmax(1)+kk1(1)*(xmax(2)-xmax(1));
  y0=ymax(1)+kk1(2)*(ymax(2)-ymax(1));
   kk2=[0.25 0.7];%
  x1=xmax(1)+kk2(1)*(xmax(2)-xmax(1));
  y1=ymax(1)+kk2(2)*(ymax(2)-ymax(1));
     kk3=[0.35 0.7];
  x2=xmax(1)+kk3(1)*(xmax(2)-xmax(1));%
  y2=ymax(1)+kk3(2)*(ymax(2)-ymax(1));
  text(x0,y0,str1,'FontSize',12,'FontWeight','bold', 'color', color_scheme_aaas(2,:),'interpreter','latex');
  text(x1,y1,str2,'FontSize',12,'FontWeight','bold', 'color', color_scheme_aaas(3,:),'interpreter','latex');
  text(x2,y2,str3,'FontSize',12,'FontWeight','bold', 'color', color_scheme_aaas(4,:),'interpreter','latex');
view(90,90)
subtitle(str4,'FontSize',25,'FontWeight','bold','interpreter','latex')
   title('(e)','position',[xmax(1)-8,ymax(1)],'FontSize',20);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%nd=1,nf=41
nexttile
for j=1:1
    numbf=linspace(1,1001,41);
    Nfobs=length(numbf);
  a_lookf=lgtt(numbf);
%Ndep=5;
Ndep=1;
%numdep=linspace(5,195,20);
numdep=[1];
Njtime1=10;
Njtime2=20;
Njtime3=20;
ksrmse=zeros(4,4);
Lx=40;
ks=K_s1(1)*ones(Nlay,Njtime2+Njtime3+Njtime1+3);
forwardmodel(K_s(:,j),K_s1,alpha,n_0,q_s,q_p);
[ks(:,1:Njtime1+1),~]=inverse_steady(ks(:,1),1,numbf,Ndep,numdep,Njtime1,Lx,K_s(:,j),q_s,q_p,alpha,n_0);
[ks(:,Njtime1+2:Njtime2+Njtime1+2),~]=inverse_amp(exp(ks(:,Njtime1+1)),Nfobs,numbf,Ndep,numdep,Njtime2,Lx,K_s(:,j),q_s,q_p,alpha,n_0);
[ks(:,Njtime2+Njtime1+3:Njtime3+Njtime2+Njtime1+3),~]=inverse_ps(exp(ks(:,Njtime2+Njtime1+2)),Nfobs,numbf,Ndep,numdep,Njtime3,Lx,K_s(:,j),q_s,q_p,alpha,n_0);
Ksestima1(:,j)=exp(ks(:,Njtime1+1));
Ksestima2(:,j)=exp(ks(:,Njtime2+Njtime1+2));
Ksestima3(:,j)=exp(ks(:,Njtime2+Njtime3+Njtime1+3));
end
cmap = get(0, 'defaultaxescolororder');
rmse0=abs(sqrt(sum((log(0.31)-log(K_s)).^2)/19)./mean(log(K_s)));
rmse1=abs(sqrt(sum((log(Ksestima1)-log(K_s)).^2)/19)./mean(log(K_s)));
rmse2=abs(sqrt(sum((log(Ksestima2)-log(K_s)).^2)/19)./mean(log(K_s)));
rmse3=abs(sqrt(sum((log(Ksestima3)-log(K_s)).^2)/19)./mean(log(K_s)));

xlab=linspace(5,195,20);
p1=plot(xlab,log(K_s),'*-','linewidth',2, 'color', color_scheme_aaas(10,:));
hold on;
p2=plot(xlab,log(Ksestima1),'*-','linewidth',2, 'color', color_scheme_aaas(2,:));
p3=plot(xlab,log(Ksestima2),'*-','linewidth',2, 'color', color_scheme_aaas(3,:));
p4=plot(xlab,log(Ksestima3),'*-','linewidth',2, 'color', color_scheme_aaas(4,:));
ylim([-2 0])
  a=get(gca);
  xmax=a.XLim;
  ymax=a.YLim;
  set(gca,'FontSize',15,'FontWeight','bold');
p5=scatter(numdep,ymax(1),'m');
    set(gca,'box','off');
set(gca,'xticklabel',[])

   str1={['$$NRMSE_1=$$' num2str(rmse1,'%4.4f') ]};
  str2={['$$NRMSE_2=$$' num2str(rmse2,'%4.4f') ]};
   str3={['$$NRMSE_3=$$' num2str(rmse3,'%4.4f') ]};
   str4={['$$n_f=$$',num2str(Nfobs),',$$n_{sp}=$$',num2str(Ndep)]};

   kk1=[0.15 0.7];%
  x0=xmax(1)+kk1(1)*(xmax(2)-xmax(1));
  y0=ymax(1)+kk1(2)*(ymax(2)-ymax(1));%
   kk2=[0.25 0.7];%
  x1=xmax(1)+kk2(1)*(xmax(2)-xmax(1));%
  y1=ymax(1)+kk2(2)*(ymax(2)-ymax(1));%
     kk3=[0.35 0.7];%
  x2=xmax(1)+kk3(1)*(xmax(2)-xmax(1));
  y2=ymax(1)+kk3(2)*(ymax(2)-ymax(1));
  text(x0,y0,str1,'FontSize',12,'FontWeight','bold', 'color', color_scheme_aaas(2,:),'interpreter','latex');
  text(x1,y1,str2,'FontSize',12,'FontWeight','bold', 'color', color_scheme_aaas(3,:),'interpreter','latex');
  text(x2,y2,str3,'FontSize',12,'FontWeight','bold', 'color',color_scheme_aaas(4,:),'interpreter','latex');
view(90,90)
subtitle(str4,'FontSize',25,'FontWeight','bold','interpreter','latex')
   title('(f)','position',[xmax(1)-8,ymax(1)],'FontSize',20);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%nd=2,nf=21
nexttile
for j=1:1
        numbf=linspace(1,501,21);
    Nfobs=length(numbf);
  a_lookf=lgtt(numbf);
%Ndep=5;
Ndep=2;
%numdep=linspace(5,195,20);
numdep=[1,101];
Njtime1=10;
Njtime2=20;
Njtime3=20;
ksrmse=zeros(4,4);
Lx=40;
ks=K_s1(1)*ones(Nlay,Njtime2+Njtime3+Njtime1+3);
forwardmodel(K_s(:,j),K_s1,alpha,n_0,q_s,q_p);
[ks(:,1:Njtime1+1),~]=inverse_steady(ks(:,1),1,numbf,Ndep,numdep,Njtime1,Lx,K_s(:,j),q_s,q_p,alpha,n_0);
[ks(:,Njtime1+2:Njtime2+Njtime1+2),~]=inverse_amp(exp(ks(:,Njtime1+1)),Nfobs,numbf,Ndep,numdep,Njtime2,Lx,K_s(:,j),q_s,q_p,alpha,n_0);
[ks(:,Njtime2+Njtime1+3:Njtime3+Njtime2+Njtime1+3),~]=inverse_ps(exp(ks(:,Njtime2+Njtime1+2)),Nfobs,numbf,Ndep,numdep,Njtime3,Lx,K_s(:,j),q_s,q_p,alpha,n_0);
Ksestima1(:,j)=exp(ks(:,Njtime1+1));
Ksestima2(:,j)=exp(ks(:,Njtime2+Njtime1+2));
Ksestima3(:,j)=exp(ks(:,Njtime2+Njtime3+Njtime1+3));
end
cmap = get(0, 'defaultaxescolororder');
rmse0=abs(sqrt(sum((log(0.31)-log(K_s)).^2)/19)./mean(log(K_s)));
rmse1=abs(sqrt(sum((log(Ksestima1)-log(K_s)).^2)/19)./mean(log(K_s)));
rmse2=abs(sqrt(sum((log(Ksestima2)-log(K_s)).^2)/19)./mean(log(K_s)));
rmse3=abs(sqrt(sum((log(Ksestima3)-log(K_s)).^2)/19)./mean(log(K_s)));
% xlabel('$$z$$(cm)','FontSize',25,'FontWeight','bold','interpreter','latex');

xlab=linspace(5,195,20);
p1=plot(xlab,log(K_s),'*-','linewidth',2, 'color', color_scheme_aaas(10,:));
hold on;
p2=plot(xlab,log(Ksestima1),'*-','linewidth',2, 'color', color_scheme_aaas(2,:));
p3=plot(xlab,log(Ksestima2),'*-','linewidth',2, 'color', color_scheme_aaas(3,:));
p4=plot(xlab,log(Ksestima3),'*-','linewidth',2, 'color', color_scheme_aaas(4,:));
ylim([-2 0])
    set(gca,'box','off');

  a=get(gca);
  xmax=a.XLim;
  ymax=a.YLim;
  set(gca,'FontSize',15,'FontWeight','bold');
p5=scatter(numdep,ymax(1),'m');

   str1={['$$NRMSE_1=$$' num2str(rmse1,'%4.4f') ]};
  str2={['$$NRMSE_2=$$' num2str(rmse2,'%4.4f') ]};
   str3={['$$NRMSE_3=$$' num2str(rmse3,'%4.4f') ]};
   str4={['$$n_f=$$',num2str(Nfobs),',$$n_{sp}=$$',num2str(Ndep)]};

   kk1=[0.15 0.7];
  x0=xmax(1)+kk1(1)*(xmax(2)-xmax(1));
  y0=ymax(1)+kk1(2)*(ymax(2)-ymax(1));
   kk2=[0.25 0.7];
  x1=xmax(1)+kk2(1)*(xmax(2)-xmax(1));
  y1=ymax(1)+kk2(2)*(ymax(2)-ymax(1));
     kk3=[0.35 0.7];%
  x2=xmax(1)+kk3(1)*(xmax(2)-xmax(1));
  y2=ymax(1)+kk3(2)*(ymax(2)-ymax(1));%
  text(x0,y0,str1,'FontSize',12,'FontWeight','bold', 'color',color_scheme_aaas(2,:),'interpreter','latex');
  text(x1,y1,str2,'FontSize',12,'FontWeight','bold', 'color', color_scheme_aaas(3,:),'interpreter','latex');
  text(x2,y2,str3,'FontSize',12,'FontWeight','bold', 'color', color_scheme_aaas(4,:),'interpreter','latex');
view(90,90)
subtitle(str4,'FontSize',25,'FontWeight','bold','interpreter','latex')
   title('(g)','position',[xmax(1)-8,ymax(1)],'FontSize',20);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%nd=2,nf=21
nexttile
for j=1:1
        numbf=linspace(501,1001,21);
    Nfobs=length(numbf);
  a_lookf=lgtt(numbf);
%Ndep=5;
Ndep=2;
%numdep=linspace(5,195,20);
numdep=[1,101];
Njtime1=10;
Njtime2=20;
Njtime3=20;
ksrmse=zeros(4,4);
Lx=40;
ks=K_s1(1)*ones(Nlay,Njtime2+Njtime3+Njtime1+3);
forwardmodel(K_s(:,j),K_s1,alpha,n_0,q_s,q_p);
[ks(:,1:Njtime1+1),~]=inverse_steady(ks(:,1),1,numbf,Ndep,numdep,Njtime1,Lx,K_s(:,j),q_s,q_p,alpha,n_0);
[ks(:,Njtime1+2:Njtime2+Njtime1+2),~]=inverse_amp(exp(ks(:,Njtime1+1)),Nfobs,numbf,Ndep,numdep,Njtime2,Lx,K_s(:,j),q_s,q_p,alpha,n_0);
[ks(:,Njtime2+Njtime1+3:Njtime3+Njtime2+Njtime1+3),~]=inverse_ps(exp(ks(:,Njtime2+Njtime1+2)),Nfobs,numbf,Ndep,numdep,Njtime3,Lx,K_s(:,j),q_s,q_p,alpha,n_0);
Ksestima1(:,j)=exp(ks(:,Njtime1+1));
Ksestima2(:,j)=exp(ks(:,Njtime2+Njtime1+2));
Ksestima3(:,j)=exp(ks(:,Njtime2+Njtime3+Njtime1+3));
end
cmap = get(0, 'defaultaxescolororder');
rmse0=abs(sqrt(sum((log(0.31)-log(K_s)).^2)/19)./mean(log(K_s)));
rmse1=abs(sqrt(sum((log(Ksestima1)-log(K_s)).^2)/19)./mean(log(K_s)));
rmse2=abs(sqrt(sum((log(Ksestima2)-log(K_s)).^2)/19)./mean(log(K_s)));
rmse3=abs(sqrt(sum((log(Ksestima3)-log(K_s)).^2)/19)./mean(log(K_s)));
  ylabel('log($$K_s$$)','FontSize',25,'FontWeight','bold','interpreter','latex');

xlab=linspace(5,195,20);
p1=plot(xlab,log(K_s),'*-','linewidth',2, 'color', color_scheme_aaas(10,:));
hold on;
p2=plot(xlab,log(Ksestima1),'*-','linewidth',2, 'color', color_scheme_aaas(2,:));
p3=plot(xlab,log(Ksestima2),'*-','linewidth',2, 'color', color_scheme_aaas(3,:));
p4=plot(xlab,log(Ksestima3),'*-','linewidth',2, 'color',color_scheme_aaas(4,:));
ylim([-2 0])

  set(gca,'FontSize',15,'FontWeight','bold');
   str1={['$$NRMSE_1=$$' num2str(rmse1,'%4.4f') ]};
  str2={['$$NRMSE_2=$$' num2str(rmse2,'%4.4f') ]};
   str3={['$$NRMSE_3=$$' num2str(rmse3,'%4.4f') ]};
   str4={['$$n_f=$$',num2str(Nfobs),',$$n_{sp}=$$',num2str(Ndep)]};
     a=get(gca);
  xmax=a.XLim;%
  ymax=a.YLim;%
  p5=scatter(numdep,ymax(1),'m');
    set(gca,'box','off');
set(gca,'xticklabel',[])

   kk1=[0.15 0.7];%
  x0=xmax(1)+kk1(1)*(xmax(2)-xmax(1));
  y0=ymax(1)+kk1(2)*(ymax(2)-ymax(1));%
   kk2=[0.25 0.7];%
  x1=xmax(1)+kk2(1)*(xmax(2)-xmax(1));
  y1=ymax(1)+kk2(2)*(ymax(2)-ymax(1));
     kk3=[0.35 0.7];
  x2=xmax(1)+kk3(1)*(xmax(2)-xmax(1));
  y2=ymax(1)+kk3(2)*(ymax(2)-ymax(1));%
  text(x0,y0,str1,'FontSize',12,'FontWeight','bold', 'color', color_scheme_aaas(2,:),'interpreter','latex');
  text(x1,y1,str2,'FontSize',12,'FontWeight','bold', 'color', color_scheme_aaas(3,:),'interpreter','latex');
  text(x2,y2,str3,'FontSize',12,'FontWeight','bold', 'color', color_scheme_aaas(4,:),'interpreter','latex');
view(90,90)
subtitle(str4,'FontSize',25,'FontWeight','bold','interpreter','latex')
   title('(h)','position',[xmax(1)-8,ymax(1)],'FontSize',20);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%nd=2,nf=41
nexttile
for j=1:1
        numbf=linspace(1,1001,41);
    Nfobs=length(numbf);
  a_lookf=lgtt(numbf);
%Ndep=5;
Ndep=2;
numdep=[1,101];
Njtime1=10;
Njtime2=20;
Njtime3=20;
ksrmse=zeros(4,4);
Lx=40;
ks=K_s1(1)*ones(Nlay,Njtime2+Njtime3+Njtime1+3);
forwardmodel(K_s(:,j),K_s1,alpha,n_0,q_s,q_p);
[ks(:,1:Njtime1+1),~]=inverse_steady(ks(:,1),1,numbf,Ndep,numdep,Njtime1,Lx,K_s(:,j),q_s,q_p,alpha,n_0);
[ks(:,Njtime1+2:Njtime2+Njtime1+2),~]=inverse_amp(exp(ks(:,Njtime1+1)),Nfobs,numbf,Ndep,numdep,Njtime2,Lx,K_s(:,j),q_s,q_p,alpha,n_0);
[ks(:,Njtime2+Njtime1+3:Njtime3+Njtime2+Njtime1+3),~]=inverse_ps(exp(ks(:,Njtime2+Njtime1+2)),Nfobs,numbf,Ndep,numdep,Njtime3,Lx,K_s(:,j),q_s,q_p,alpha,n_0);
Ksestima1(:,j)=exp(ks(:,Njtime1+1));
Ksestima2(:,j)=exp(ks(:,Njtime2+Njtime1+2));
Ksestima3(:,j)=exp(ks(:,Njtime2+Njtime3+Njtime1+3));
end
cmap = get(0, 'defaultaxescolororder');
rmse0=abs(sqrt(sum((log(0.31)-log(K_s)).^2)/19)./mean(log(K_s)));
rmse1=abs(sqrt(sum((log(Ksestima1)-log(K_s)).^2)/19)./mean(log(K_s)));
rmse2=abs(sqrt(sum((log(Ksestima2)-log(K_s)).^2)/19)./mean(log(K_s)));
rmse3=abs(sqrt(sum((log(Ksestima3)-log(K_s)).^2)/19)./mean(log(K_s)));

xlab=linspace(5,195,20);
p1=plot(xlab,log(K_s),'*-','linewidth',2, 'color', color_scheme_aaas(10,:));
hold on;
p2=plot(xlab,log(Ksestima1),'*-','linewidth',2, 'color', color_scheme_aaas(2,:));
p3=plot(xlab,log(Ksestima2),'*-','linewidth',2, 'color', color_scheme_aaas(3,:));
p4=plot(xlab,log(Ksestima3),'*-','linewidth',2, 'color', color_scheme_aaas(4,:));
  set(gca,'FontSize',15,'FontWeight','bold');
   str1={['$$NRMSE_1=$$' num2str(rmse1,'%4.4f') ]};
  str2={['$$NRMSE_2=$$' num2str(rmse2,'%4.4f') ]};
   str3={['$$NRMSE_3=$$' num2str(rmse3,'%4.4f') ]};
   str4={['$$n_f=$$',num2str(Nfobs),',$$n_{sp}=$$',num2str(Ndep)]};
   ylim([-2 0])
     a=get(gca);
  xmax=a.XLim;%
  ymax=a.YLim;%
  p5=scatter(numdep,ymax(1),'m');
    set(gca,'box','off');
set(gca,'xticklabel',[])

   kk1=[0.15 0.7];%
  x0=xmax(1)+kk1(1)*(xmax(2)-xmax(1));%
  y0=ymax(1)+kk1(2)*(ymax(2)-ymax(1));
   kk2=[0.25 0.7];%
  x1=xmax(1)+kk2(1)*(xmax(2)-xmax(1));%
  y1=ymax(1)+kk2(2)*(ymax(2)-ymax(1));%
     kk3=[0.35 0.7];%
  x2=xmax(1)+kk3(1)*(xmax(2)-xmax(1));%
  y2=ymax(1)+kk3(2)*(ymax(2)-ymax(1));%
  text(x0,y0,str1,'FontSize',12,'FontWeight','bold', 'color', color_scheme_aaas(2,:),'interpreter','latex');
  text(x1,y1,str2,'FontSize',12,'FontWeight','bold', 'color', color_scheme_aaas(3,:),'interpreter','latex');
  text(x2,y2,str3,'FontSize',12,'FontWeight','bold', 'color', color_scheme_aaas(4,:),'interpreter','latex');
view(90,90)
subtitle(str4,'FontSize',25,'FontWeight','bold','interpreter','latex')
   title('(i)','position',[xmax(1)-8,ymax(1)],'FontSize',20);
h1=legend([p1,p2,p3,p4,p5],'True value','Assimilate $$\theta_{c}$$','Assimilate $$\theta_{c}$$+$$Amp$$','Assimilate $$\theta_{c}$$+$$Amp$$+$$PS$$','Observation','interpreter','latex');
h1.Location='southeast';
h1.FontSize = 25;
set(h1,'Box','off')
 tt.TileSpacing = 'compact';
 tt.Padding = 'compact';
xlabel(tt,'log($$K_s$$)','FontSize',25,'FontWeight','bold','interpreter','latex')
ylabel(tt,'$$z$$(cm)','FontSize',25,'FontWeight','bold','interpreter','latex')
h1.Layout.Tile = 'east';

set(gcf, 'PaperSize', [50*1.2 30*1.2])
saveas(gcf,'Fig3.pdf')