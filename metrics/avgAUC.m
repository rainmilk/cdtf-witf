function [avgScore, aucscore] = avgAUC(X, positiveIdx, negativeIdx)
% avgAUC
%   X: Data
%   positiveIdx: indecies for postive instance
%   negativeIdx: indecies for negative instances

nUser = size(positiveIdx,1);
aucscore = zeros(nUser,1);
parfor n = 1:nUser
    pIdx = positiveIdx(n,:);
    Xn = X(n,:);
    if any(pIdx)
        nIdx = negativeIdx(n,:);
        pScore = Xn(pIdx);
        nScore = Xn(nIdx);
        aucscore(n) = auc(pScore, nScore);
    end
end
avgScore = mean(aucscore(aucscore>0));
end
