function meanXk = meandata(Xk, implicit)
if implicit
    nXk = numel(Xk);
else
    nXk = nnz(Xk);
end
meanXk = sum(nonzeros(Xk))/nXk;
end