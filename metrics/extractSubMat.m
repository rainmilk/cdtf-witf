function [m] = extractSubMat( rank, rankrange )
    [~,idx] = sort(rank, 'descend');
    start = max(1, ceil(rankrange(1)));
    endt = min(length(idx), ceil(rankrange(2)));
    m = idx(start:endt);
end

