function results = m_hist_Isot(R,ml,mmax,ret)
%% Obtención de m-entornos direccionales.
%
% R    = Número de observaciones.
% ml   = longitud de la m-historia en cada dirección.
% mmax = es una cota superior de control para no incluir en la m-historia
% vecinos extremadamente lejanos. Busca los ml mas próximos en cada una de
% las direcciones entre los mmax mas próximos.
% ret = tipo de reticula (regular=cuadrada o irregular=randn).
%
% SALIDA: mhi = m-historia direccional del cuadrante i-ésimo.
%   Se dan los m entornos de aquellas localizaciones que tienen m-entorno
%   completo en todas las direcciones.
%
% Utiliza: COOR ; m_hist0
%% Selección del tipo de reticula.
if ret==1,
C=COOR(R,1);C=C.C;
elseif ret==0,
    C=randn(2,R);
end
%% Especificaciones
mh=m_hist0(C(:,1),C(:,2),mmax);
C1=C(:,1);
C2=C(:,2);
C1mh=C1(mh);
C2mh=C2(mh);
for i=2:mmax
    D1(:,i-1)=C1mh(:,i)-C1mh(:,1);
    D2(:,i-1)=C2mh(:,i)-C2mh(:,1);
end

%% mh1: mh del primer cuadrante
for j=1:R
for i=1:mmax-1
    if D1(j,i)>0 && D2(j,i)>=0, D3(j,i)=1; 
    elseif D1(j,i)<=0 && D2(j,i)>0, D3(j,i)=2; 
    elseif D1(j,i)<0 && D2(j,i)<=0, D3(j,i)=3;
    elseif D1(j,i)>=0 && D2(j,i)<0, D3(j,i)=4;
    end
end
end
mh1=zeros(R,mmax);
mh2=zeros(R,mmax);
mh3=zeros(R,mmax);
mh4=zeros(R,mmax);
for i=1:R
    k1=0;k2=0;k3=0;k4=0;
    for j=1:mmax-1
        if D3(i,j)==1, k1=k1+1;mh1(i,k1)=mh(i,j+1);
        elseif D3(i,j)==2, k2=k2+1;mh2(i,k2)=mh(i,j+1);
        elseif D3(i,j)==3, k3=k3+1;mh3(i,k3)=mh(i,j+1);
        elseif D3(i,j)==4, k4=k4+1;mh4(i,k4)=mh(i,j+1);
        end
    end
end
%ml=4;
mh1=mh1(:,1:ml);
mh2=mh2(:,1:ml);
mh3=mh3(:,1:ml);
mh4=mh4(:,1:ml);
%% Selección de observaciones simbolizables
Front=sum((mh1==0)*1,2)+sum((mh2==0)*1,2)+sum((mh3==0)*1,2)+sum((mh4==0)*1,2);
ff=sum((Front==0));
mhf1=zeros(ff,ml+1);mhf2=zeros(ff,ml+1);mhf3=zeros(ff,ml+1);mhf4=zeros(ff,ml+1);
kk=0;
for i=1:R
    if Front(i)==0,kk=kk+1;mhf1(kk,:)=[i mh1(i,:)];mhf2(kk,:)=[i mh2(i,:)];mhf3(kk,:)=[i mh3(i,:)];mhf4(kk,:)=[i mh4(i,:)];end
end
results.mh1=mhf1;
results.mh2=mhf2;
results.mh3=mhf3;
results.mh4=mhf4;