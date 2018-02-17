function [Xk, Wi] = imputenoise(Xk, Wk, implicit, rnParam)
[N,Mk] = size(Wk);
cs = sum(logical(Wk),2) < rnParam.Tau;
ncs = nnz(cs);
nimpute = min(rnParam.Size, N);
if nimpute*ncs > 0.2*numel(Wk)
    Wi = false(size(Wk));
else
    Wi = logical(spalloc(N, Mk, nimpute*ncs));
end

mu = rnParam.Mean;
if isempty(mu)
    mu = meandata(Xk,implicit);
end
Std = rnParam.Std;

parfor i=1:N
    if cs(i)
        idx = false(1,Mk);
        idx(randperm(Mk,nimpute)) = true;
        idx(logical(Wk(i,:))) = false;
        Wi(i,:) = idx;
        Xgi = Xk(i,:);
        Xgi(idx) = mu + Std*randn(1,nnz(idx));
        Xk(i,:) = Xgi;
    end
end
end

