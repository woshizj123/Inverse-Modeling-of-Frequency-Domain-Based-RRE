%% Townley number Convert to frequency f
function [ff]=q_switch(Ks,alpha,n_0,lgtt)

v=((Ks)./(n_0));

% lgf=linspace(-3,2,1001);
% f=10.^(lgf);
% for i=1:1001
% omega(i)=2*pi()*f(i);
% tt(i)=sqrt((sqrt(1+16*(omega(i)^2)/((k_p*alpha(1))^2))+1)/2);
% end
% aa=log10(tt-1);
% lgtt=linspace(-3,1,1001);
aa1=10.^(lgtt)+1;
ff=sqrt(((2*aa1.^2-1).^2-1)*v^2*alpha^2/16)/2/pi();

