%function [Emcs,Emcr,Ehmms,Ehmmr]=Conditional_DHMM_KernelEntropy(mAtractor)
%
% This function estimates en conditional hidden Markov entropy of a embedded
% attractor. The estimation is based on a discrete hidden Markov Model with 
% non-parametric estimation of the distributions on each state of the model
% using the Parzen window method and a clustering using a method based on 
% Principle of Relevant Information.



function [Salida,mAtractor_cuantizado]=Conditional_DHMM_KernelCorrEntropy2(mAtractor,iLong,iAlpha,rParam,flags,mCentroids)
tiempo = tic;

if nargin < 1, error('mAtrator is required'); end
if nargin < 2, iLong = size(mAtractor,1); end
if nargin < 3, iAlpha = 2; end
if nargin < 4, rParam=0.2; end
if nargin < 5, flags = 'MHR'; end

mAtractor2 = mAtractor(1:iLong,:);
iDim = size(mAtractor,2);
%----------------------------------------------------------------------
%Silverman's sigma
S = mAtractor2;
deviation = std(S);
sigmax = mean(deviation);
sigma = sigmax*((4/(size(S,1)*(2*iDim+1)))^(1/(iDim+4)));
%------------------------------------------------------------------
if nargin < 6
    %mCentroids=clustering_PRI(mAtractor2);
    [mCentroids] = SCMS(mAtractor2,sigma);
end
%-------------------------------------------------------------------
%------------- Eliminar centroides repetidos -----------------------
mCentroids2 = mCentroids;
for i=1:size(mCentroids,1)
    vCentroid = mCentroids2(i,:);
    mDist = mCentroids2 - ones(size(mCentroids2,1),1)*vCentroid;
    mDist = sqrt(sum(mDist.^2,2));
    ind = mDist <= 1e-5; % umbral arbitrario % 0.3; % umbral arbitrario
    ind(i)=0;
    mCentroids2(ind,:)=[];
    if size(mCentroids2,1) <= i
        break;
    end
end

mCentroids = PrunePC(mCentroids2,mAtractor2);
%-------------------------------------------------------------------
mAtractor_cuantizado = Quantization(mCentroids,mAtractor);
%--------------------------------------------------------------------------
%                Calcular matriz de transición
[mTransmat,~,mCentroids,bFlag] = dmm_transmat(mAtractor_cuantizado,mCentroids);
if bFlag == 1
    mAtractor_cuantizado=Quantization(mCentroids,mAtractor);
end
if isempty(mAtractor_cuantizado)
    %pause();
    error('Attractor is empty');
end
%--------------------------------------------------------------------------
if any( flags=='M' )
    pos = (mTransmat == 0);
    mTransmat2 = mTransmat + pos;
    %--------------------------------------------------------------------------
    tem = mTransmat2.*log2(mTransmat2);
    Salida.Emcs = -sum(tem(:))/size(mTransmat,1);
    %--------------------------------------------------------------------------
    tem = mTransmat.^iAlpha;
    tem2 = sum(tem,2);
    pos = (tem2 == 0);
    tem2 = log2(tem2 + pos)/(1-iAlpha);
    Salida.Emcr = sum(tem2)/size(mTransmat,1);
end
%--------------------------------------------------------------------------
if any( flags=='R' )
    [~,~,Salida.Ecorrs,Salida.Ecorrr,Salida.MRDs,Salida.MRDr]=MarkovRecurrenceEntropy(mAtractor,mAtractor_cuantizado,mTransmat);
end
if any( flags=='H' )
    [Eps,Epr,SigSilv,iEfectState]=ParzenEntropy(mAtractor,mAtractor_cuantizado,iAlpha,rParam);
    %--------------------------------------------------------------------------
    Salida.Ehmms = Eps/iEfectState; 
    Salida.Ehmmr = Epr/iEfectState;
    Salida.Ehmms2 = Eps/10;
    Salida.Ehmmr2 = Epr/10;
    Salida.Ehmmsig = SigSilv/iEfectState;
end


toc(tiempo)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Julián David Arias Londoño %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Universidad de Antioquia, Colombia %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Septiembre de 2012 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
