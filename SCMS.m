function [mYSal] = SCMS(mX,h)
%function [PC] = SCMS(X,h) estimates the Principal Curve using the method
% Subspace Constrained Mean Shift, see:
% Ozertem, Umut; Erdogmus, Deniz, "Locally Defined Principal Curves and
% Surfaces", Journal of Machine Learning Research 12 (2011) 241-274
%
%
plotFlag = 0;
if nargin < 2, h = 0.01; end %Sería mejor usar Sigma de silverman


%Initizalization of method's variables
kEpsilon = 0.005;
kd = 1;
%--------------------------------------------------------------------------
%Initizalization of Algorithm's variables
[iN,iD] = size(mX);
c = (2*pi)^(-iD/2);
mYSal = zeros(1,iD);
%rBandSq = 5*h;
rBandSq = 0.61;
%rBandSq =  -0.2*h/4.2 + 0.3;
j=0;
iMaxIter = 100;
beenVisitedFlag = zeros(1,iN,'uint8');              %track if a points been seen already
PuntosNoVisitados = iN;
initPtInds = find(beenVisitedFlag == 0);
%--------------------------------------------------------------------------
while PuntosNoVisitados ~= 0
    %tem = randperm(iN);
    vY = mX(initPtInds(1),:);
    fFlag = 1;
    iIter = 0;
    myMembers = [];
    j = j + 1;
    while iIter<iMaxIter && fFlag == 1
        %j = j+1;
        iIter = iIter + 1;
        mY = repmat(vY,iN,1);
        %mYSal(j,:) = vY;
        %mDK = Dkernel_mat(mX,vY,h,'gauss');
        mK = kernel_mat(mX,vY,h,'gauss');
        inInds = find(mK > rBandSq);
        if ~isempty(inInds)
            myMembers   = [myMembers inInds];
            beenVisitedFlag(myMembers) = 1;
            vNum = mX(inInds,:).*repmat(mK(inInds)',1,iD);
            iDen = sum(mK(inInds));
            vM = sum(vNum,1)/iDen - vY; %m
            %----------------------------------------------------------------------
            vNum2 = mX(inInds,:).*repmat(mK(inInds)',1,iD);
            iDen2 = sum(mK(inInds));
            vM2 = sum(vNum2,1)/iDen2;
            %----------------------------------------------------------------------
            vG = sum((mY(inInds,:)-mX(inInds,:)).*repmat(mK(inInds)',1,iD),1)/(iN*(h^2));
            rP = sum(mK)/iN;
            mH = ((((mY(inInds,:)-mX(inInds,:))*2.*repmat(mK(inInds)',1,iD))'./(h^2))*(mY(inInds,:)-mX(inInds,:)))*(c/(iN*h^(2+iD))) - sum(mK(inInds)).*eye(iD)*c/(iN*h^(2+iD));
            mISigma = -(1/rP)*mH + (1/(rP^2))*(vG')*vG;
            [V,D] = eig(mISigma);
            vEigVal = diag(D);
            [~,vInd] = sort(vEigVal,1,'descend');
            PC = V(:,vInd(1:iD-kd))';
            vY2 = (PC')*PC*(vM') + vY';
            if plotFlag
                figure(1200),clf,hold on
                if iD == 2
                    plot(mX(:,1),mX(:,2),'.')
                    plot(mX(myMembers,1),mX(myMembers,2),'ys')
                    plot(vM2(1),vM2(2),'go')
                    plot(vY(1),vY(2),'rd')
                    plot(vY2(1),vY2(2),'kd')
                    pause(0.3)
                end
            end
            if norm(vY2'-vY) < kEpsilon
                fFlag = 0;
            else
                vY = vY2';
            end
        else
            fFlag = 0;
        end
    end
    mYSal(j,:) = vY;
    initPtInds = find(beenVisitedFlag == 0);           %we can initialize with any of the points not yet visited
    PuntosNoVisitados = length(initPtInds);                   %number of active points in set
end

