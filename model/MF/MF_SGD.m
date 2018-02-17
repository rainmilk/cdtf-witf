function [U V] = MF_SGD(Y,m,n,U,V,lambda,lr,decrate,maxiter,tol)
%% [U V] = MF_SGD(Y,m,n,U,V,lambda,lr,decrate,maxiter,tol)
% Inputs:
%   Y: observation matrix
%   [m, n] = index pairs of Obseravations
%   U: latent factor matrix, e.g. initialized randomly U = abs(randn(size(Y,1),50))
%   V: latent factor matrix, e.g. initialized randomly V = abs(randn(size(Y,2),50))
%   lamda: regularization parameters
%   lr: learning rate
%   decrate: decrease rate for lr
%   maxiter: the number of maximum iterations
%   tol: tolerance for checking convergences
% Outputs:
%   U: learned latent factor matrix
%   V: learned latent factor matrix

%% default values
if nargin<6
    lambda=0.02;
end
if nargin<7
    lr=0.005;
end
if nargin<8
    decrate=1;
end
if nargin<9
    maxiter=1000;
end
if nargin<10
    tol=1e-5;
end

%%
res = 0;
nval = length(m);
globalt = tic;
for iter=1:maxiter
    oldres = res;
    res = 0;
    tic;
    p = randperm(nval);
    m = m(p);
    n = n(p);
    for k=1:nval
        i=m(k);
        j=n(k);
        oldVal = U(i,:)*V(j,:)';
        err = Y(i,j) - oldVal;
        U(i,:) = U(i,:) + lr * (err * V(j,:) - lambda * U(i,:));
        V(j,:) = V(j,:) + lr * (err * U(i,:) - lambda * V(j,:));
        newVal = U(i,:)*V(j,:)';
        res = res + 0.5*(Y(i,j) - newVal)^2;
    end
    lr = lr*decrate;
    tolval = abs(res - oldres);
    fprintf('\nIteration %d: %g secs, Residual: %g, Tolerance: %g.', iter, toc, res, tolval);
    if tolval < tol
        break;
    end
end

fprintf('\nFinished in %d iterations, %g secs.\n', iter, toc(globalt));