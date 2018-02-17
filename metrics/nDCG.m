function [ score ] = nDCG( rel, irel )
% NDCG
%   rel: rank for related instances
%   irel: rand for irrelated instances
n=length(rel);
DCG = rel(1);
iDCG = irel(1);
if n>1
    iDCG = iDCG + sum(irel(2:n)./log2(2:n));
    DCG = DCG + sum(rel(2:n)./log2(2:n));
end
score=DCG/iDCG;
end

