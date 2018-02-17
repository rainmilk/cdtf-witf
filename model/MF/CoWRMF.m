function [ U, V ] = CoWRMF( X, W, lambda_U, lambda_V, U, V, muU, muV, alpha, varargin )
% nirank Neighborhood based Implicit Feedback Ranking
%   Detailed explanation goes here
%   L = ai*sum_ij ||[Wi(i,j)] * [Xi(i,j) - U(i,:)*Vi(j,:)']||^2
%           + lambda_U*||U - muU||^2  + cplx*lambda_U*||U||^2
%           + lambda_Vi*||Vi - muVi||^2 + cplx*lambda_Vi*||Vi||^2
%
%   U(i,:) = ai * ( Xi(i,:)*diag(Wi(i,:))*V1 + lambda_U*muU(i,:))*inv(Vi'*diag(W(i,:))*Vi+lambda_U*I)
%
%   V(j,:) = (X(:,j)'*diag(W(:,j))*U + + lambda_V*muV(i,:))*inv(U'*diag(W(:,j))*U+lamba_V*I)



%% Set parameters
nM = length(X);

params = inputParser;
params.addParamValue('MaxIt',200, @isscalar); % Max number of iteration
params.addParamValue('AbsTol', 1e-10, @isscalar);  % Fit error to exit
params.addParamValue('RelTol', 1e-6, @isscalar);  % Fit error to exit
params.addParamValue('implicit', false, @isscalar);  % If implicit feedback
params.addParamValue('epsilon', 0.01*ones(1,nM), @isnumeric);  % epsilon for IF
params.addParamValue('CplxCtrl', 1, @isscalar);  % Control complexity, cplx
params.parse(varargin{:});

MaxIt   = params.Results.MaxIt;
AbsTol = params.Results.AbsTol;
RelTol = params.Results.RelTol;
implicit = params.Results.implicit;
epsilon = params.Results.epsilon;
cplx = params.Results.CplxCtrl;

[nU, R] = size(U);
X_W = cell(1,nM);
nV = zeros(1,nM);
for m=1:nM
    nV(m) = size(V{m},1);
    X_W{m} = X{m} .* W{m};
    if implicit, X_W{m} = X_W{m} + epsilon*X{m}; end
end

%% Set regularization parameters
if nargin < 3 || isempty(lambda_U)
    lambda_U = 0;
end

if nargin < 4 || isempty(lambda_V)
    lambda_V = zeros(1,nM);
end

%% Set prior means
if nargin < 7 || isempty(muU), muU = spalloc(nU,R,0); end
if nargin < 8 || isempty(muV)
    muV = cell(1,nM);
    for m=1:nM
        muV{m} = spalloc(nV(m),R,0);
    end
end;

lambda_muU = lambda_U * muU;
lambda_U = cplx * lambda_U * speye(R);

lambda_muV = cell(1,nM);
for m=1:nM
    lambda_muV{m} = lambda_V(m) * muV{m};
end

%% Iterations
UtU = 0;
VtV = 0;
oldfit = inf;
tid = tic;
for it=1:MaxIt
    
    %% Update UVtV
    Uold = U;
    
    U = zeros(size(U));
    for m=1:nM
        Wm = W{m};
        XWm = X_W{m};
        Vm = V{m};
        alpham = alpha(m);
        if implicit, VtV = epsilon * (Vm'*Vm); end
        parfor i=1:nU
            Wi = Wm(i,:);
            maski = logical(Wi);
            XWi = XWm(i,:);
            lambda_muUi = lambda_muU(i,:);
            Vi = Vm(maski,:);
            Wi = Wi(maski);
            A = XWi(maski) * Vi + lambda_muUi;
            B = Vi' .* repmat(Wi,R,1) * Vi + lambda_U;
            if implicit, B = B + VtV; end
            U(i,:) =  U(i,:) + alpham * (A / B);
        end
    end
    
    if implicit, UtU = epsilon * (U'*U); end
    %% Update V
    Vold = V;
    for m=1:nM
        Wm = W{m};
        XWm = X_W{m};
        Vm = V{m};
        lambda_muVm = lambda_muV{m};
        lambda_Vm = cplx * lambda_V(m) * speye(R);
        parfor j=1:nV(m)
            Wj = Wm(:,j);
            maskj = logical(Wj);
            XWj =  XWm(:,j);
            lambda_muVj = lambda_muVm(j,:);
            Uj = U(maskj,:);
            Wj = Wj(maskj);
            A = XWj(maskj)' * Uj + lambda_muVj;
            B = (Uj .* repmat(Wj,1,R))' * Uj + lambda_Vm;
            if implicit, B = B + UtU; end
            Vm(j,:) = A / B;
        end
        V{m} = Vm;
    end
    
    fit = norm(U-Uold, 'fro');
    for m=1:nM
        fit = fit + norm(V{m}-Vold{m}, 'fro');
    end
    tol = abs(fit-oldfit)/oldfit;
    
    fprintf('\nFactor Variation (AbsTol): %g\tIter: %g\tRelTol: %g\tExpire: %g secs',fit,it,tol,toc(tid));
    
    if tol <= RelTol || fit <= AbsTol
        break;
    end
    
    oldfit = fit;
end
end
