function results = QMap(Y,mh,Simb,info)
% PURPOSE: computes QMap statistic for Maps Comparations
% ---------------------------------------------------
%  USAGE: results = GS_se(Y,mh,Simb,plot_cond)
%  where: Y  = T*n vector of qualitative variables entered as 1,2,...,K
%         mh = n*m matrix with the m-history of each spatial observation
%         symb = (k^m)*m matrix containing symbols
%         info = an (optional) structure variable with input options
%         info.plot  = 'yes' / ['no'] flag for plotting symbols histogram
%         info.print = 'yes' / ['no'] print the result
% ---------------------------------------------------
%  RETURNS: results.gs - the value of the statistic
%           results.S  - vector for producing histogram with frequency of symbols
%  --------------------------------------------------
%  SEE ALSO: m_hist0 ; m_hist_no
% ---------------------------------------------------
% REFERENCES: 
% ---------------------------------------------------

% written by:
% Fernando A. López Hernández, 30/06/2009
% Facultad de C.C. de la Empresa
% Dpto. Metodos Cuantitativos e Informáticos
% Universidad Politécnica de Cartagena
% E-mail: Fernando.Lopez@upct.es

% Modificacion del test GS para m-entornos no solapados
[T N]=size(Y);
[R m]=size(mh);
Ymax=max(max(Y));
%[n1 m]=size(mh);
[nusi n2]=size(Simb);
lns=zeros(1,nusi);
S=zeros(1,R);
%N=length(Y);
%if length(Y)~=n1, error('Error: mh and Y are different lengths'), end 
%if m~=n2, error('Error: my and symb are different lengths'), end 
if nargin<4
    info.plot='no';
    info.print='no';
else
    info.plot=lower(info.plot);
end
%% Simbolizar la serie de datos.
for t=1:T
    A=Y(t,:);
    Z(:,:,t)=A(mh);
end
for t=1:T
for i=1:R
    for s=1:nusi
        %Z(i,:,t)
        if isequal(Simb(s,:),Z(i,:,t))==1, S(t,i)=s; end
    end
end
results.S=S;
end
for t=1:T
nsk(t,:)=hist(S(t,:),1:nusi);
end
%nsk

if strcmp(info.plot,'yes')
    for f=1:T
    figure(f);hist(S(f,:)',1:nusi);
    end
end

fnk=nsk/(R*T);
fnT=mean(nsk)/(R*T);
for h=1:T
for i=1:nusi
    if nsk(h,i)==0, lnsk(h,i)=0; else lnsk(h,i)=log(fnT(i)/fnk(h,i)); end
end
end
QMap=-2*sum(sum(nsk.*lnsk));
results.QMap=QMap;
results.Rns=R/nusi;
results.QMap2=[QMap 1-chi2cdf(QMap,(T-1)*nusi) chi2inv(0.95,(T-1)*nusi) (T-1)*nusi];

if strcmp(info.print,'yes')
fprintf(1,'\n')
fprintf(1,'  Sample size               (N) = %9.0f\n',N)
fprintf(1,'  Cortes Temporales         (T) = %9.0f\n',T)
fprintf(1,'  Observacines Simbolizadas (R) = %9.0f\n',R)
fprintf(1,'  Número de modalidades     (k) = %9.0f\n',Ymax)
fprintf(1,'  Longitud m-historia       (m) = %9.0f\n',m)
fprintf(1,'  Grado de solapamiento     (r) = %9.0f\n',m-(N-m)/(R-1))
fprintf(1,'  Número de Símbolos        (n) = %9.0f\n',nusi)
fprintf(1,'  Ratio R/n                     = %9.2f\n',R/nusi)
%fprintf(1,'  Probabilidades =                ');fprintf(1,'%9.4f ',pm');fprintf(1,'\n');
fprintf(1,'\n')
fprintf(1,'Test QMapS for Maps comparation qualitative data\n');
fprintf(1,'        Test                      Value     DF    p-value \n');
fprintf(1,'        QM  =          %9.2f',QMap);fprintf(1,'    %2.0f',(T-1)*nusi);fprintf(1,' %9.4f \n',1-chi2cdf(QMap,(T-1)*nusi));
end

