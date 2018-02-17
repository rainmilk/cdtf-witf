function [ result ] = meanPrecision2( X, positiveIdx, negativeIdx,atK)
%AVGRECALL Summary of this function goes here
%   Detailed explanation goes here
result = zeros(length(atK), 1);
nUser = size(X,1);
nTest = 0;
for n = 1:nUser
    pIdx = positiveIdx(n,:);
    if nnz(pIdx)>0
        nTest = nTest + 1;
        nIdx = negativeIdx(n,:);
        pScore = X(n, pIdx);
        nScore = X(n, nIdx);
        result = result + precisionAtK(pScore, nScore, atK);
    end
end
if nTest > 0
    result = result / nTest;
end
end

