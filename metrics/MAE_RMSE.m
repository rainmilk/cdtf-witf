function [MAE, RMSE, nR, allerrs] = MAE_RMSE(U,V,R,rank,rankrange)

Rs = nonzeros(R);
maxR = max(Rs);
minR = min(Rs);
clear Rs;

nGrp = 1;

emptyrank = false;
if nargin > 4 && ~isempty(rankrange)
    nGrp = length(rankrange);
    if isempty(rank)
        emptyrank = true;
    end
end

outputerr = false;
allerrs = cell(1, nGrp);
if nargout > 3
    outputerr = true;
end

MAE = zeros(1, nGrp);
RMSE = zeros(1, nGrp);
nR = zeros(1, nGrp);

Rg = R;
Vg = V;

for g = 1:nGrp
    if nGrp > 1
        if emptyrank
            subidx = rankrange{g};
        else
            subidx = extractSubMat(rank, rankrange{g});
        end
        Rg = R(:,subidx);
        Vg = V(:,subidx);
    end
    
    M = size(Rg,2);
    MAEg = 0;
    RMSEg = 0;
    
    errs = [];
    if outputerr
        errs = cell(1, M);
    end
    
    parfor j=1:M
        Rj = Rg(:,j);
        mask = logical(Rj);
        pred = Vg(:,j)'*U(:, mask);
        pred(pred>maxR) = maxR;
        pred(pred<minR) = minR;
        err = Rj(mask)' - pred;
        if outputerr
            errs{j} = err;
        end
        
        MAEg = MAEg + sum(abs(err));
        RMSEg = RMSEg + err*err';
    end
    nR(g) = nnz(Rg);
    MAE(g) = MAEg/nR(g);
    RMSE(g) = sqrt(RMSEg/nR(g));
    
    if outputerr
       allerrs{g} = errs;
    end
end
end