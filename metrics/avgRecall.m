function [ result ] = avgRecall( X, positiveIdx, negativeIdx, atK)
% AVGRECALL
%   X: Data
%   positiveIdx: indecies for postive instance
%   negativeIdx: indecies for negative instances
%   atK: recall@K

result = zeros(1, length(atK));
nUser = size(X,1);
nTest = 0;
parfor n = 1:nUser
    pIdx = positiveIdx(n,:);
    Xn = X(n,:);
    if nnz(pIdx)>0
        nTest = nTest + 1;
        nIdx = negativeIdx(n,:);
        pScore = Xn(pIdx);
        nScore = Xn(nIdx);
        result = result + recallAtK(pScore, nScore, atK);
    end
end
if nTest > 0
    result = result / nTest;
end
end

