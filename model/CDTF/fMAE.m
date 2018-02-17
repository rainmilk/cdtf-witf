function [ mae, A,H,C,P ] = fMAE( wX, cache, R, lambda, X, testIdx, testMat, conv, ...
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
    mae = cache(weights);
    fprintf('fMAE Id: %g, find MAE: %g', id, mae);
else
    if ~exist('params','var')
        params.dummy = [];
    end
    [A,H,C,P,fit]=wparafac2(X,R,Constraint,Opt,lambda,wX,params);
    mae = predict;
    
    it = 0;
    while tol>conv
        it = it + 1;
        params.initFac = {A,H,C,P};
        [A,H,C,P,fit]=wparafac2(X,R,Constraint,Opt,lambda,wX,params);
        mae_i = predict;
        tol = (mae - mae_i)/mae;
        if mae_i < mae
            mae = mae_i;
        end
        fprintf('fMAE Id: %g iteration %d MAE: %g, Tol: %g\n', id, it, mae_i, tol);
    end
    cache(weights) = mae;
end
fprintf('Leave fMAE Id: %g\n', id);

    function [MAE,RMSE] = predict
        UT = bsxfun(@times, A, C(1,:));
        VT = P{1}*H;
        [MAE,RMSE] = MAE_RMSE(UT', VT', testMat);
    end
end
