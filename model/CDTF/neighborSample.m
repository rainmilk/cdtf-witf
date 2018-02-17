function [ Xk ] = neighborSample( X, neighbor, nExample )
[N,M] = size(X);
Xk = zeros(size(X));
nimpute = min(nExample,M);
parfor i=1:N
    nbData = X(neighbor(i,:),:);
    nbIdx = logical(nbData);
    nbIdx(:, logical(X(i,:))) = false;
    idx = find(nbIdx);
    len = length(idx);
    idx = idx(randperm(len, min(nimpute,len)));
    x = Xk(i,:);
    [~,j] = ind2sub(size(nbData),idx);
    x(j) = nbData(idx);
    Xk(i,:) = x;
end
end

