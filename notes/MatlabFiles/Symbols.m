function results = Symbols(k,m)
%% Obtención de los simbolos: Datos Cualitativos
maple('with(combinat)');
cs=[];
for j=1:k
for i=1:m
cs=strcat(cs,',',num2str(j));
end
end
cs=cs(2:length(cs));
eval(['s=maple(''permute([',cs,'],',num2str(m),')'');'])
s2=str2num(s);
s2=reshape(s2',m,k^m)';
results=s2;
