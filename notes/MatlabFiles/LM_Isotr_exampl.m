%% LM_Isotr
clear
N=900;
R=N^.5;
XI=randn(N,1);
a0=2;
b0=3;
WI=W_Direcc(R);
W1=WI.W1;
W2=WI.W2;
W3=WI.W3;
W4=WI.W4;
W=WI.W;

WS=W2+W4;
WE=W1+W3;

A=inv(eye(N,N)-0.0*normw(W)-0.9*normw(WS));
e=randn(N,1);
y=A*(a0+b0.*XI+e);
x=[ones(N,1) XI];
figure(2),surfc(reshape(y,R,R))
LM_Isotr(y,x,R,2)

load Simbolos
info.plot='no';
info.print='yes';
ml=2;
mmax=20;
m=m_hist_Isot(N,ml,mmax,1);
Isot(y',m,Simb0,info);
