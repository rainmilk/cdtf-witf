function [ score ] = auc(pScore, nScore, atK )
% AUC
%   positiveScore: scores for postive instances
%   negativeScore: scores for negative instances
M=length(pScore);
N=length(nScore);
[~,IX] = sort([pScore, nScore]);
ranking(IX) = 1:(M+N);
pRank = ranking(1:M);
score=(sum(pRank)-(M*(M+1)/2))/(M*N);

if nargin > 2
    score=repmat(score,1,length(atK));
end

end
