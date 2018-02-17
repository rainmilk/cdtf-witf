function [ recall, hits, T ] = recallAtK( pScore, nScore, atK )
% recallAtK
%   pScore: scores for positive instances
%   nScore: scores for negative instances
%   atK:    top K

T = length(pScore);
N = length(nScore);
R = T + N;
[~,IX] = sort([pScore, nScore],'descend');
rankvec = [ones(1,T), zeros(1,N)];
rankvec = cumsum( rankvec(IX) );

topK = min(atK,R);
hits = rankvec(topK);
recall = hits/T;

