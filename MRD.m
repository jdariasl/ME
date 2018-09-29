function [Shannon,Renyi2]=MRD(vTemporalDist)

flag = 1;
i = 0;
iAlpha = 2;
while flag
    i = i +1;
    Elemento = vTemporalDist(1);
    RepElemento = find(vTemporalDist == Elemento);
    vNElementos(i) = length(RepElemento);
    vTemporalDist(RepElemento)=[];
    if isempty(vTemporalDist)
        flag = 0;
    end
end

iNElementos = sum(vNElementos);
vNElementos = vNElementos/iNElementos;
tem = vNElementos.*log2(vNElementos);
Shannon = -sum(tem);
Renyi2 = (1/(1-iAlpha))*log2(sum(vNElementos.^iAlpha));

    