function [ rmse, A,H,C,P ] = fRMSE( wX, cache, R, lambda, X, testIdx, trueVals, conv, ...
    Constraint, Opt, params)
% fMAE 
%   wX: weights (not including the 1st slice)
%   R:  the number of factors
%   lambda: regularization parameter
%   X:  Data

id = rand;
fprintf('Entering fMAE Id: %g\n', id);
tol = realmax;
wX = [1,wX];
weights = mat2str(wX);
fprintf('Weights: %s\n', weights);
if cache.isKey(weights)
    rmse = cache(weights);
    fprintf('fRMSE Id: %g, find MAE: %g\n', id, rmse);
else
    if ~exist('params','var')
        params.dummy = [];
    end
    [A,H,C,P,fit]=wparafac2(X,R,Constraint,Opt,lambda,wX,params);
    rmse = predict;
    
    it = 0;
    while tol>conv
        it = it + 1;
        params.initFac = {A,H,C,P};
        [A,H,C,P,fit]=wparafac2(X,R,Constraint,Opt,lambda,wX,params);
        rmse_i = predict;
        tol = (rmse - rmse_i)/rmse;
        if rmse_i < rmse
            rmse = rmse_i;
        end
        fprintf('fRMSE Id: %g iteration %d MAE: %g, Tol: %g\n', id, it, rmse_i, tol);
    end
    cache(weights) = rmse;
end
fprintf('Leave fMAE Id: %g\n', id);

    function rmse = predict
        predvalues = A*diag(sparse(C(1,:)))*(P{1}*H)';
        predvalues = predvalues(testIdx);
        predvalues( predvalues > 5 ) = 5;
        predvalues( predvalues < 1 ) = 1;
        rmse = sqrt(mean(( trueVals - predvalues).^2));
    end
end
