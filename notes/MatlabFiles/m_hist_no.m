function results = m_hist_no(xc,yc,m,s)
% PROPOSITO: Obtiene las m-historias eligiendo los m-1 vecinos mas próximos
% a cada localizacion con índice de solapamiento s.
% ---------------------------------------------
% USAGE: 
% where: xc = coordenadas x.
%        yc = coordenadas y.
%        m  = tamaño de la historia.
%        s  = indice de solapamiento.
% ---------------------------------------------
%
% Escrita en fecha 07/01/09 por:
%
% Fernando A. López Hernández
% Facultad de C.C. de la Empresa
% Dpto. Métodos Cuantitativos e Informáticos
% Universidad Politécnica de Cartagena
% E-mail: fernando.lopez@upct.es
% 
%-----------------------------------------------
n=length(xc);
nn=find_neighbors(xc,yc,n-1);
nnlist=zeros(n,m-1);
ns=fix((n-m)/(m-s))+1;
list=zeros(1,ns);
list(1)=1;
nnlist(1,:)=nn(1,1:m-1);
listanegra=[1 nnlist(1,1:m-(s+1))];
t=1;
for v=2:ns
list(v)=nnlist(t,m-s);%list(v)=nnlist(t,m-1);
h=list(v);
a=nn(h,:);
k=0;
for i=1:n-1
    if sum(listanegra==a(i))==0,k=k+1;nnlist(h,k)=a(i);t=h;end
    if k==m-1, break, end
end
listanegra=[listanegra h nnlist(h,1:m-(s+1))];
end
nnlist=([(1:n)' nnlist]);
v=0;
for i=1:n
    if nnlist(i,2)~=0, v=v+1; nnlist1(v,:)=nnlist(i,:);
    end
end
results=nnlist1;


