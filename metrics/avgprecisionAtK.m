function [ avp ] = avgprecisionAtK( pScore, nScore, atK )
T = length(pScore);
N = length(nScore);
R = T + N;

if nargin<3
    atK = R;
end

[~,IX] = sort([pScore, nScore],'descend');
rankvec = [ones(1,T), zeros(1,N)];
rankvec = rankvec(IX);


topK = min(atK,R);
avp = (cumsum(rankvec).*rankvec) ./ (1:R);
avp = cumsum(avp);
avp = avp(topK) ./ min(topK,T);