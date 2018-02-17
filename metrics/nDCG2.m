function [ ndcg ] = nDCG2( pScore, nScore, atK )
% NDCG
%   rel: rank for related instances
%   irel: rand for irrelated instances

T = length(pScore);
N = length(nScore);

den = log2( 2:(T+N+1) );

rankvec = [ones(1,T), zeros(1,N)];
IDCG = cumsum( (2.^rankvec - 1) ./  den);

[~,IX] = sort([pScore, nScore],'descend');
rankvec = rankvec(IX);
DCG = cumsum( (2.^rankvec - 1) ./ den );

topK = min( atK, T+N );


ndcg = DCG(topK) ./ IDCG(topK);

end

