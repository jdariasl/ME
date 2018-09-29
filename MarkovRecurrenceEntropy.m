function [Emcs,Emcr,Ecorrs,Ecorrr,MRDs,MRDr]=MarkovRecurrenceEntropy(mAtractor,mAtractor_cuantizado,mTransmat)
%iDim = size(mAtractor,2);
iAlpha = 2;
%--------------------------------------------------------------------------
pos = (mTransmat == 0);
mTransmat2 = mTransmat + pos;
%--------------------------------------------------------------------------
tem = mTransmat2.*log2(mTransmat2);
%tem2 = eModelo.vPrior*ones(1,Nstate);
%tem = tem.*tem2; % Como se entrena con una �nica secuencia, Pi = 1,0,0,0
Emcs = -sum(sum(tem))/size(mTransmat,1);
%--------------------------------------------------------------------------
tem = mTransmat.^iAlpha;
tem2 = sum(tem,2);
pos = (tem2 == 0);
tem2 = log2(tem2 + pos)/(1-iAlpha);
Emcr = sum(tem2)/size(mTransmat,1);
%--------------------------------------------------------------------------
Nstates = max(mAtractor_cuantizado);
vEps = zeros(1,Nstates);
vEpr = vEps;
vMRDs = vEpr;
vMRDr = vEpr;
for i=1:Nstates
    vInd = find(mAtractor_cuantizado==i);
    if ~isempty(vInd)
        vTemporalDist = abs(diff(vInd(end:-1:1)));
        if length(vTemporalDist) > 1
            %Estima entropias de shannon y renyi de orden para la recuencia del
            %estado i
            [vMRDs(i),vMRDr(i)]=MRD(vTemporalDist);
%             %--------------------------------------------------------------
%             %Encuentra subsecuencias dentro del estado i
%             mSequences = FindSecuences(vTemporalDist,vInd);
%             %Calcula la Matriz Kernel basada en Correlacion
%             mKernelCorrMatrix = KernelCorrMatrix(mAtractor,mSequences);
%             iN = size(mKernelCorrMatrix,1);
%             %----------------------------------------------------------
%             tem = sum(mKernelCorrMatrix,2)/iN;
%             temR = (iN*tem).^(iAlpha-1);
%             tem = log2(tem);
%             vEps(i) = -sum(tem)/iN;
%             vEpr(i) = (1/(1-iAlpha))*log2(sum(temR)/(iN^iAlpha));
            %--------------------------------------------------------------
            vEps(i) = 0;
            vEpr(i) = 0;
        else
            vMRDs(i) = 0;
            vMRDr(i) = 0;
            vEps(i) = 0;
            vEpr(i) = 0;
        end
    else
        sText = strcat('The estate_',num2str(i),'_does not exist!');
        disp(sText);
    end
end
Ecorrs = sum(vEps)/Nstates;
Ecorrr = sum(vEpr)/Nstates;
MRDs = sum(vMRDs)/Nstates;
MRDr = sum(vMRDr)/Nstates;
        
    
    
    