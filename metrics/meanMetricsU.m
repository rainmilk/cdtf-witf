function [ varargout ] = meanMetricsU( metrics, U, V, posIdx, negIdx, atK, isTrIdx, rank, rankrange)
% meanMetrics
%   U, V: factor matrices for user and item
%   train, test: training and test set
%   atK: recall@K
nGrp = 1;

if nargin < 7, isTrIdx = false; end
emptyrank = false;
if nargin > 7
    nGrp = length(rankrange);
    if isempty(rank)
        emptyrank = true;
    end
end

nMetric = length(metrics);

K = max(length(atK), 1);

varargout = cell(nMetric, 1);
for i = 1:nMetric
    varargout{i} = zeros(nGrp, K);
end

posIdxR = posIdx;
negIdxR = negIdx;
VR = V;
nTester = size(posIdxR, 1);
for g = 1:nGrp
    if nGrp > 1
        if emptyrank
            subidx = rankrange{g};
        else
            subidx = extractSubMat(rank, rankrange{g});
        end
        posIdxR = posIdx(:,subidx);
        negIdxR = negIdx(:,subidx);
        VR = V(:,subidx);
    end
    
    output = zeros(nMetric, K);
    nTest = 0;
    parfor n = 1:nTester
        pIdx = posIdxR(n,:);                            % positive idx
        nIdx = negIdxR(n,:);                            % negative idx
        if isTrIdx, nIdx = ~( pIdx | nIdx ); end;       % if negIdx is traning idx
        if any(pIdx) && any(nIdx)
            nTest = nTest + 1;
            pScore = U(:,n)' * VR(:,pIdx);
            if isTrIdx
                nScore = UR(:,n)' * V;
                nScore = nScore(nIdx);
            else
                nScore = UR(:,n)' * V(:,nIdx);
            end
            result = zeros(nMetric, K);
            for i=1:nMetric
                metric = metrics{i};
                result(i,:) = metric(pScore, nScore, atK);
            end
            output = output + result;
        end
    end
    output = output ./ nTest;
    for i=1:nMetric
        varargout{i}(g,:) = output(i,:);
    end
end
end

