function results = Qsemi_variog(Y,C,di,Simb,ret)
%% PROPOSITO: Obtención de semivariogramas empíricos direccionales para variables cualitativas.
% Semivariograma en dos direcciones. N-S y E-O.
% también se pueden hacer 4 direcciones.
% Se recomienda obtener solo hasta h=mmd la mitad de la distancia máxima.

% ENTRADAS:
% Y = Vector de datos R2x1.
% C = Coordenadas de cada punto. Vector R2x2
% di = vector 1x2 que marca las dos direcciones de los variogramas en grados
%     entre [0,360] con tolerancia +/- 45º
% ret = selecciona tipo reticula (reg=1/irreg=2) 
%       selecc mejor los valores de h.
% SALIDA: Figura con semivariogramas
%
% VER TAMBIEN: semi_variog, COOR, wt, wr, make_neighborsw
%
% Escrita en fecha 14/8/09 por:
%
% Fernando A. López Hernández
% Facultad de C.C. de la Empresa
% Dpto. Metodos Cuantitativos e Informaticos
% Universidad Politécnica de Cartagena
% e-mail: Fernando.Lopez@upct.es

%% 
R2=length(C);

%% Matriz de distancias entre puntos
d=zeros(R2,R2);
x=C(:,1);
y=C(:,2);
for i=1:R2
   for j=1:R2
      d(i,j)=sqrt((x(i)-x(j))^2+(y(i)-y(j))^2);
   end
end
%% TT(i,j) indica el ángulo que forma la obs. j con repescto a la obs. i. 
[T,R]=cart2pol((kron(x',ones(R2,1))-kron(x,ones(1,R2))),(kron(y',ones(R2,1))-kron(y,ones(1,R2))));
TT=(T/(2*pi))*360;
TT=(TT>=0).*TT+(TT<0).*(360+TT);

%% semivariograma omnidireccional y direccionales
% Elige los h.
if ret==1, % retucula regular.
    mmd=(R2^.5)/2;
    h=(0:mmd);
elseif ret==2, %reticula irregular.
mmd=R2^.5;
mm=max(max(d))/2;
h=[0:mm/mmd:mm];
mmd=10;
end

%% Matrices direccionales
% % Para cuatro direcciones [0,45,90,135] 
% da(1).eq=(TT>=0).*(TT<45)+(TT>=180).*(TT<180+45);
% da(2).eq=(TT>=45).*(TT<45+45)+(TT>=45+180).*(TT<45+180+45);
% da(3).eq=(TT>=90).*(TT<90+45)+(TT>=90+180).*(TT<90+180+45);
% da(4).eq=(TT>=135).*(TT<135+45)+(TT>=135+180).*(TT<135+180+45);

% Para dos direcciones generalmente [0,90] y [90,180] 
tol=45; % la tolerancia de la direccion
d1=di(1);
d2=di(2);

% da(1).eq=(TT>=d1-tol).*(TT<d1+tol)+(TT>=d1-tol+180).*(TT<d1+tol+180);
% da(2).eq=(TT>=d2-tol).*(TT<d2+tol)+(TT>=d2-tol+180).*(TT<d2+tol+180);

if d1==45,
da(1).eq=(TT>=0).*(TT<90)+(TT>=180).*(TT<270);
da(2).eq=(TT>=90).*(TT<180)+(TT>=270).*(TT<360);
elseif d1==0
da(1).eq=(TT>=360-45).*(TT<45)+(TT>=90+45).*(TT<45+180);
da(2).eq=(TT>=45).*(TT<90+45)+(TT>=45+180).*(TT<270+45);
end


[nusi n2]=size(Simb);

%% Q Omnivariograma.

for i=1:mmd
    
    dd=((((d>h(i)).*d).*((d<=h(i+1)).*d))>0)*1;% dd(i,j)=1 sii distancia(obs i;obs j) esta entre [h(i),h(i+1)]
    ff=(dd.*ones(R2,R2)).*kron(Y,ones(1,R2))'; 
    S=[kron(Y,ones(R2,1)),reshape(ff',R2*R2,1)];
    SS=sortrows(S,-2);
    nd=sum(sum(dd.*ones(R2,R2)));
    SS=SS(1:nd,:);

    % Asigna simbolos
Sa=zeros(nd,1);
for k=1:nusi
    Sa=Sa+ismember(SS,Simb(k,:),'rows')*k;
end

    nsk=hist(Sa,1:nusi);
    NSK(i,:)=nsk;
    for j=1:nusi
    if nsk(j)==0, lns(j)=0; else lns(j)=log(nsk(j)/nd); end
    end
    Ohm(i)=-sum((nsk/nd).*lns);
    clear Sa
    clear lns
    clear nsk

end
clear dd
clear ff
clear S
clear SS
clear Sa
clear lns
clear nsk

%% Q Semivariograma dirección 1.
for i=1:mmd
    dd=((((d>h(i)).*d).*((d<=h(i+1)).*d))>0)*1;% dd(i,j)=1 sii distancia(obs i;obs j) esta entre [h(i),h(i+1)]
    ff=(dd.*da(1).eq).*kron(Y,ones(1,R2))'; 
    S=[kron(Y,ones(R2,1)),reshape(ff',R2*R2,1)];
    SS=sortrows(S,-2);
    nd=sum(sum(dd.*da(1).eq));
    SS=SS(1:nd,:);

    % Asigna simbolos
Sa=zeros(nd,1);
for k=1:nusi
    Sa=Sa+ismember(SS,Simb(k,:),'rows')*k;
end

    nsk=hist(Sa,1:nusi);
    NSK1(i,:)=nsk;
    for j=1:nusi
    if nsk(j)==0, lns(j)=0; else lns(j)=log(nsk(j)/nd); end
    end
    hm(i,1)=-sum((nsk/nd).*lns);ns1=nsk;
    clear Sa
    clear lns
    clear nsk

end
clear dd
clear ff
clear S
clear SS
clear Sa
clear lns
clear nsk
%% Q Semivariograma dirección 2.

for i=1:mmd
    dd=((((d>h(i)).*d).*((d<=h(i+1)).*d))>0)*1;% dd(i,j)=1 sii distancia(obs i;obs j) esta entre [h(i),h(i+1)]
    ff=(dd.*da(2).eq).*kron(Y,ones(1,R2))'; 
    S=[kron(Y,ones(R2,1)),reshape(ff',R2*R2,1)];
    SS=sortrows(S,-2);
    nd=sum(sum(dd.*da(2).eq));
    SS=SS(1:nd,:);
    
% Asigna símbolos.     
Sa=zeros(nd,1);
for k=1:nusi
    Sa=Sa+ismember(SS,Simb(k,:),'rows')*k;
end

    nsk=hist(Sa,1:nusi);
    NSK2(i,:)=nsk;
    for j=1:nusi
    if nsk(j)==0, lns(j)=0; else lns(j)=log(nsk(j)/nd); end
    end
    hm(i,2)=-sum((nsk/nd).*lns);
    clear Sa
    clear lns
    clear nsk

end
%%

figure(1),subplot(1,2,1),plot(Ohm,'-'),
xlabel('h')
ylabel('\gamma(h)')
title('Omnivariograma')

if max(Y)==2, minx=min(min(hm))-0.1; elseif max(Y)==3, minx=min(min(hm))-0.1; end
axis([0 mmd minx log(nusi)+0.1])
subplot(1,2,2),plot(hm,'-')
title('Semivariogramas Direccionales')
xlabel('h')
ylabel('\gamma(h)')
axis([0 mmd minx log(nusi)+0.1])
if ret==1,
figure(2),pcolor(reshape(Y,R2^.5,R2^.5))
else
    a=['+b';'*g';'^r';'+y'];
    for i=1:max(Y)
    figure(2),plot(C(:,1).*(Y==i) ,C(:,2).*(Y==i),a(i,:))
    hold on
    end
    hold off
end

hm
results.NSK=NSK;
results.NSK1=NSK1;
results.NSK2=NSK2;
return
