function [mTransmat,mTransmatOrg,mCentroids,bFlag] = dmm_transmat(datos_cuant,mCentroids)
iNcO = size(mCentroids,1);
iTrans = length(datos_cuant);
mTransmat = zeros(max(datos_cuant),max(datos_cuant));
for i=2:iTrans
    mTransmat(datos_cuant(i-1),datos_cuant(i)) = mTransmat(datos_cuant(i-1),datos_cuant(i)) + 1;
end
mTransmatOrg = mTransmat;
%mTransmat = mTransmat + 1; % Regla de laplace
vSuma = sum(mTransmat,2);
flag=min(vSuma);
while flag == 0
    ind = vSuma ~= 0;
    mTransmat = mTransmat(ind,ind);
    mCentroids = mCentroids(ind,:);
    vSuma = sum(mTransmat,2);
    flag=min(vSuma);
end
mTransmat = mTransmat./(vSuma*ones(1,size(mTransmat,1)));
if iNcO ~= size(mCentroids,1);
    bFlag = 1;
else
    bFlag = 0;
end
%--------------------------------------------------------------------------