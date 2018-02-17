function [ aucscore, A,H,C,P] = fAUC( wX, cache, R, alpha, beta, lambda, X,  positiveIdx, negativeIdx, conv, ...
    Constraint, Opt, params)
% fAUC 
%   wX: weights (not including the 1st slice)
%   R:  the number of factors
%   lambda: regularization parameter
%   X:  Data

id = rand;
fprintf('Entering fAUC Id: %g\n', id);
tol = realmax;
wX = [1,wX];
weights = mat2str(wX);
fprintf('Weights: %s\n', weights);
if cache.isKey(weights)
    aucscore = cache(weights);
    fprintf('fAUC Id: %g, find AUC: %g', id, aucscore);
else
    if ~exist('params','var')
        params.dummy = [];
    end
    
    X = if_input( X, wX, alpha, beta);
    
    [A,H,C,P,fit]=wparafac2(X,R,Constraint,Opt,lambda,wX,params);
    aucscore = avgAUC(A*diag(sparse(C(1,:)))*(P{1}*H)',positiveIdx,negativeIdx);
    fprintf('fAUC Id: %g iteration 0 AUC: %g\n', id, aucscore);
    
    it = 0;
    while tol>conv
        it = it + 1;
        params.initFac = {A,H,C,P};
        [A,H,C,P,fit]=wparafac2(X,R,Constraint,Opt,lambda,wX,params);
        aucscore_i = avgAUC(A*diag(sparse(C(1,:)))*(P{1}*H)',positiveIdx,negativeIdx);
        tol = (aucscore_i - aucscore)/aucscore;
        if aucscore_i > aucscore
            aucscore = aucscore_i;
        end
        fprintf('fAUC Id: %g iteration %d AUC: %g, Tol: %g\n', id, it, aucscore_i, tol);
    end
    cache(weights) = aucscore;
end
fprintf('Leave fAUC Id: %g\n', id);
aucscore = -aucscore;
end
