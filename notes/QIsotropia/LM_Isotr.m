function results = LM_Isotr(y,x,Wf,C,cond)
% PROPOSITO: Obtención del test LM de isotropía.
% Conceptualmente es el mismo que el test LM del papel de Dijon 
%            modelo SAR con expresión funcional para la dep espacial.
%  La toolbox se llama LM_casetti.m en Millave/Paris08.
%  La Manga 8 de Agosto 2009. 
%
% ---------------------------------------------
% USAGE: 
% where: y = variable dependiente.
%        x = vector de variables independientes.
%        xc = vector de coordenadas (latitud).
%        yc = vector de coordenadas (longitud).
%        w = matriz de contactos. w en principio estandarizada.
%     cond = toma el valor i para un polinomio de grado i (i=1,2).
% ---------------------------------------------
%
% Escrita en fecha 01/12/08 por:
%
% Fernando A. López Hernández
% Facultad de C.C. de la Empresa
% Dpto. Metodos Cuantitativos e Informaticos
% Universidad Politécnica de Cartagena
% E-mail: fernando.lopez@upct.es
% 
%-----------------------------------------------

[n sx]=size(x);%n=length(y);
I=speye(n);
if cond==1,
    np=3;
else
    np=5;
end
%w=normw(w);

%% Selección de las matrices de Isoptropia W_Direcc ó W_Isotr
WI=W_Direcc4(Wf,C);

W(1).eq=normw(WI.W);
w=normw(WI.W);
W(2).eq=normw(WI.W1);
W(3).eq=normw(WI.W2);
W(4).eq=normw(WI.W3);
W(5).eq=normw(WI.W4);

% Definicion opciones para el proceso optimizacion:
options0 = optimset('fminbnd');

% minimizacion -Lik
% Obtencion del minimo bajo la hip nula.
mh0=fminbnd('f_exph0',-1,1,options0,y,x,w);

% Evalua modelo bajo hipotesis nula.
rho=mh0;
B=I-rho*w;
Beta=inv(x'*x)*(x'*B*y);
u=B*y-x*Beta;
sigma2=u'*u/n;
results.rho=rho;
results.Beta=Beta';
results.sigma2=sigma2;
% 
BI=inv(B);

%% Evaluar la matriz de varianzas y covarianzas bajo hipotesis alternativa.

e11=x'*x; g=length(e11);
e12=x'*W(1).eq'*BI*x*Beta;
for i=2:np
e12=horzcat(e12,x'*W(i).eq'*BI*x*Beta);
end
e13=zeros(g,1);
e21=e12';
e22=zeros(np,np);
A=BI*x*Beta;
for i=1:np
    for j=1:np
    e22(i,j)=A'*W(i).eq'*W(j).eq*A + sigma2*(trace(BI*W(i).eq*BI*W(j).eq)+trace(BI*W(i).eq'*BI'*W(j).eq));
    end
end
e23=trace(BI*W(1).eq);
for i=2:np
e23=vertcat(e23,trace(BI*W(i).eq));
end
e31=e13';
e32=e23';
e33=n/(2*sigma2);
MI=(1/sigma2)*[e11 e12 e13;e21 e22 e23;e31 e32 e33];

%% Gradiente
e=(B*y-x*Beta);
gr2=(e'*W(1).eq*y-sigma2*trace(BI*W(1).eq));
for i=2:np
gr2=vertcat(gr2,(e'*W(i).eq*y-sigma2*trace(BI*W(i).eq)));
end
gr=(1/sigma2)*[zeros(1,sx) gr2' 0];
LM=gr*inv(MI)*gr';
results.LM=LM;