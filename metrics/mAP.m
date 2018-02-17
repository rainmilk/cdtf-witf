function [score,scores] = mAP(actual, prediction, k)
% MAP: Calculates the average precision at k
%   score = mAP(actual, prediction, k)
%
%   actual is a cell array of vectors
%   prediction is a cell array of vectors
%   k is an integer
%

scores = zeros(length(prediction),1);

for i=1:length(prediction)
    if nargin<3
        scores(i) = avgprecision(actual{i}, prediction{i});
    else
        scores(i) = avgprecision(actual{i}, prediction{i}, k);
    end
end

score = mean(scores);