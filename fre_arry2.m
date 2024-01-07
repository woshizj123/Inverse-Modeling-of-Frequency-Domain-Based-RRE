%% This is a function of choosing the frequency position
function [N,L_arry]=fre_arry2(N)

i=1;
L_arry(1)=1;
while  i<N
    L_arry(i+1)=L_arry(i)+25;
    i=i+1;
end
