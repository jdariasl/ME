function [Eps,Epr,SigSilv,iEfectState]=ParzenEntropy(mAtractor,mAtractor_cuantizado,iAlpha,rParam)

nStates = max(mAtractor_cuantizado);
iDim = size(mAtractor,2);
mEps = zeros(1,nStates);
mEpr = mEps;
mSigSilv = mEpr;
iEfectState = 0;
for i=1:nStates
    
    ind = (mAtractor_cuantizado == i);
    mAtractor2 = mAtractor(ind,:);
    N = size(mAtractor2,1);
    if N > 0
        iEfectState = iEfectState +1;
        S = mAtractor2;
        %----------------------------------------------------------------------
        %Estimación de sigma de Silverman 
        if nargin < 4
            deviation = std(S);
            sigmax = mean(deviation);
            sigma = sigmax*((4/(size(S,1)*(2*iDim+1)))^(1/(iDim+4)));
        else
            sigma = rParam;
        end
        %-----------------------------------------------------------rma-----------
        norma = 2;
        if norma == 2
            K = kernel_mat(S,S,sigma,'gauss');
        elseif norma == 3 %norma infinito
            K = kernel_mat2(S,S,sigma,'gauss');
        end
        tem2 = isinf(K);
        K(tem2) = 1;
        tem = sum(K,2)/(N);
        temR = (N*tem).^(iAlpha-1);
        tem = log2(tem);
        mEps(i) = -sum(tem)/N;
        mEpr(i) = (1/(1-iAlpha))*log2(sum(temR)/(N^iAlpha));
        mSigSilv(i) = sigma;
    end
end
Eps = sum(mEps);
Epr = sum(mEpr);
SigSilv = sum(mSigSilv);
disp(strcat('Number of states in the model =',num2str(iEfectState)));    