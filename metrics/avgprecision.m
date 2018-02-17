function score = avgprecision(actual, prediction, k)
%AVERAGEPRECISIONATK   Calculates the average precision at k
%   score = avgprecision(actual, prediction, k)
%
%   actual is a vector
%   prediction is a vector
%   k is an integer
%

if nargin<3 || k>length(prediction)
    k=length(prediction);
end

score = 0;
num_hits = 0;
num_rel = min(length(actual), k);
for i=1:k
    if any(actual==prediction(i))
        num_hits = num_hits + 1;
        score = score + num_hits / i;
        if num_hits>= num_rel
            break;
        end
    end
end

score = score / num_rel;