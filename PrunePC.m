 function mCentroids = PrunePC(mCentroids2,mAtractor2)
    
    while size(mAtractor2,1) > 100 && size(mCentroids2,1) > 0.1*size(mAtractor2,1)
        %disp('Entró a podado');
        [N,D] = size(mCentroids2);
        tem = zeros(N,D);
        tem(1:end-1,:) = mCentroids2(2:end,:);
        tem(end,:) = mCentroids2(1,:);
        Dist = sqrt(sum((mCentroids2 - tem).^2,2));
        [~,ind] = min(Dist);
        if ind < N
            medio = (mCentroids2(ind,:) + mCentroids2(ind+1,:))/2;
            mCentroids2(ind,:) = medio;
            mCentroids2(ind + 1,:) = [];
        else
            medio = (mCentroids2(ind,:) + mCentroids2(1,:))/2;
            mCentroids2(ind,:) = medio;
            mCentroids2(1,:) = [];
        end
    end
    mCentroids = mCentroids2;

        
        