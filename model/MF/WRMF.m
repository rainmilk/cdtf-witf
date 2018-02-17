function [ U, V ] = WRMF( X, W, lambda_U, lambda_V, U, V, muU, muV, varargin )
% nirank Neighborhood based Implicit Feedback Ranking
%   Detailed explanation goes here
%   L = sum_ij ||[W(i,j)] * [X(i,j) - U(i,:)*V(j,:)']||^2
%              +lambda_U*||U - muU||^2 + lambda_V*||V - muV||^2
%              +cplx*lambda_U*||U||^2 + cplx*lambda_V*||V||^2
%
%   U(i,:) = ( X(i,:)*diag(W(i,:))*V + lambda_U*muU(i,:)) * inv(V'*diag(W(i,:))*V+lambda_U*I)
%
%   V(j,:) = (X(:,j)*diag(W(:,j))*U + + lambda_V*muV(i,:)) * inv(U'*diag(W(:,j))*U+lamba_V*I)



%% Set parameters
params = inputParser;
params.addParamValue('MaxIt',200, @isscalar); % Max number of iteration
params.addParamValue('AbsTol', 1e-10, @isscalar);  % Fit error to exit
params.addParamValue('RelTol', 1e-6, @isscalar);  % Fit error to exit
params.addParamValue('implicit', false, @isscalar);  % If implicit feedback
params.addParamValue('epsilon', 0.01, @isscalar);  % epsilon for IF
params.addParamValue('CplxCtrl', 1, @isscalar);  % Control complexity, cplx
params.parse(varargin{:});

MaxIt   = params.Results.MaxIt;
AbsTol = params.Results.AbsTol;
RelTol = params.Results.RelTol;
implicit = params.Results.implicit;
epsilon = params.Results.epsilon;
cplx = params.Results.CplxCtrl;

[nU, R] = size(U);
nV = size(V,1);

mask = logical(W);

lU = any(mask,2);
nlU = nnz(~lU);
lV = any(mask,1)';
nlV = nnz(~lV);

%% Set regularization parameters
if nargin < 3 || isempty(lambda_U)
    lambda_U = 0;
end

if nargin < 4 || isempty(lambda_V)
    lambda_V = 0;
end

%% Set prior means
if nargin < 7 || isempty(muU), muU = spalloc(nU,R,0); end
if nargin < 8 || isempty(muV), muV = spalloc(nV,R,0); end

lambda_muU = lambda_U * muU;
lambda_muV = lambda_V * muV;

lambda_U = cplx * lambda_U * speye(R);
lambda_V = cplx * lambda_V * speye(R);

X_W = X .* W;
if implicit, X_W = X_W + epsilon * X; end

%% Iterations
UtU = 0;
VtV = 0;
oldfit = inf;
tid = tic;
for it=1:MaxIt
    
    %% Update UVtV
    Uold = U;
    if implicit, VtV = epsilon * (V'*V); end
    parfor i=1:nU
        if lU(i)
            maski = mask(i,:);
            XWi =  X_W(i,:);
            XWi = XWi(maski);
            Wi = W(i,:);
            Wi = Wi(maski);
            lambda_muUi = lambda_muU(i,:);
            Vi = V(maski,:);
            A = XWi * Vi + lambda_muUi;
            B = repmat(Wi,R,1) .* Vi' * Vi + lambda_U;
            if implicit, B = B + VtV; end
            U(i,:) =  A / B;
        end
    end
    U(~lU,:) = repmat(mean(U(lU,:)),nlU,1);
    
    %% Update V
    Vold = V;
    if implicit, UtU = epsilon * (U'*U); end
    parfor j=1:nV
        if lV(j)
            maskj = mask(:,j);
            XWj =  X_W(:,j);
            XWj = XWj(maskj);
            Wj = W(:,j);
            Wj = Wj(maskj);
            lambda_muVj = lambda_muV(j,:);
            Uj = U(maskj,:);
            A = XWj' * Uj + lambda_muVj;
            B = (repmat(Wj,1,R) .* Uj)' * Uj + lambda_V;
            if implicit, B = B + UtU; end
            V(j,:) = A / B;          
        end
    end
    V(~lV,:) = repmat(mean(V(lV,:)),nlV,1);
    
    fit = norm(U-Uold, 'fro') + norm(V-Vold, 'fro');
    tol = abs(fit-oldfit)/oldfit;
    
    fprintf('\nFactor Variation (AbsTol): %g\tIter: %g\tRelTol: %g\tExpire: %g secs',fit,it,tol,toc(tid));
    
    if tol <= RelTol || fit <= AbsTol
        break;
    end
    
    oldfit = fit;
end

end

