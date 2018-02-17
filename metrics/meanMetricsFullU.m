function [ varargout ] = meanMetricsFullU( metrics, X, posIdx, negIdx, atK, rank, rankrange)
% AVGRECALLEX
%   U, V: factor matrices for user and item
%   train, test: training and test set
%   atK: recall@K
nGrp = 1;

emptyrank = true;
if nargin > 6
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
XR = X;
nTester = size(posIdxR, 1);
for g = 1:nGrp
    if nargin > 5
        if emptyrank
            subidx = rankrange{g};
        else
            subidx = extractSubMat(rank, rankrange{g});
        end
        posIdxR = posIdx(:,subidx);
        negIdxR = negIdx(:,subidx);
        XR = X(:,subidx);
    end
    
    output = zeros(nMetric, K);
    nTest = 0;
    parfor n = 1:nTester
        pIdx = posIdxR(n,:);                            % positive idx
        nIdx = negIdxR(n,:);                            % negative idx
        Xn = XR(n,:);
        if any(pIdx) && any(nIdx)
            nTest = nTest + 1;
            pScore = Xn(pIdx);
            nScore = Xn(nIdx);
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
