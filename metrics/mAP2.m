function [score,scores] = mAP2(X, actInd, testInd, k)
% MAP: Calculates the average precision at k
%   score = mAP2(X, actInd, preInd, k)
%
% 

scores = zeros(1,size(X,1));
for i=1:length(scores)
    actual = find(actInd(i,:));
    if ~isempty(actual)
        rec = X(i,:);
        [~,IX] = sort(rec,'descend');
        testIdx = testInd(i,:);
        predict = IX(testIdx(IX));
        if nargin > 3
            scores(i) = avgprecision(actual,predict,k);
        else
            scores(i) = avgprecision(actual,predict);
        end
    end
end

score = mean(scores(scores>0));

end

