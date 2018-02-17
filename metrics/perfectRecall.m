function [ avgRecall ] = perfectRecall( posIdx, negIdx, atK, rank, rankrange)
% AVGRECALLEX
%   U, V: factor matrices for user and item
%   train, test: training and test set
%   atK: recall@K

if nargin > 3
    items = extractSubMat(rank, rankrange);
    posIdx = posIdx(:,items);
    negIdx = negIdx(:,items);
end

%recall = zeros(length(atK), 1);
nUser = size(negIdx, 1);
hits = zeros(length(atK), nUser);
T = zeros(1, nUser);
parfor n = 1:nUser
    nP = nnz(posIdx(n,:));
    nN = nnz(negIdx(n,:));
    if nP>0 && nN>0
        [~, hits(:,n), T(n)] = recallAtK(ones(1,nP), zeros(1,nN), atK);
    end
end

avgRecall = mean( bsxfun( @times, hits(:, T>0), 1./T(T>0) ), 2) ;
end

