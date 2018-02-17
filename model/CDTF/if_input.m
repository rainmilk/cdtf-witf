function [ X ] = if_input( X, w, alpha, beta )
%IF_INPUT
%  Construct the reparameterized matrices for implicit feedback model
%   X: rating matrices
%   w: weight for each slice
%   alpha: scale factor
%   beta: basic offset

if nargin < 4
    beta = 50;
end

parfor k=1:length(X)
    Y = X{k};
    X{k} = w(k) * (alpha*Y + beta*logical(Y));
end
end

