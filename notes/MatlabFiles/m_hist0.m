function results = m_hist0(xc,yc,m)
% PROPOSITO: Obtiene las m-historias eligiendo los m-1 vecinos mas pr�ximos
% a cada localizacion y ordenandolas en el sentido de las agujas del reloj
% empezando a las 12.
% ---------------------------------------------
% USAGE: 
% where: xc = coordenadas x.
%        yc = coordenadas y.
%        m  = tama�o de la historia.
% ---------------------------------------------
%
% Escrita en fecha 17/01/08 por:
%
% Fernando A. L�pez Hern�ndez
% Facultad de C.C. de la Empresa
% Dpto. Metodos Cuantitativos e Inform�ticos
% Universidad Polit�cnica de Cartagena
% E-mail: Fernando.Lopez@upct.es
% 
%-----------------------------------------------

n=length(xc);
nnlist = zeros(n,m); %aqu� pondre las n m-historias

for i=1:n;
    xi = xc(i);
    yi = yc(i);
dist = (xc - xi*ones(n,1)).^2 + (yc - yi*ones(n,1)).^2;
% Introduce una precisi�n maxima de 2 decimales
dist=fix(dist*100)/100;
for j=1:n
    [THETA(j),RHO(j)] = cart2pol( -xc(j)+xi,-yc(j)+yi);
end
    THETAT=THETA';
    DC=[dist THETAT];
    [xds xind] = sortrows(DC,[1 -2]);
nnlist(i,1:m) = xind(2:m+1,1)';
end;
aa=1:n;
sal=[aa' nnlist];
results=sal;
