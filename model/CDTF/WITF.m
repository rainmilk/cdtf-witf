function [U,V,C,P,fit]=WITF(X, R, W, wX, implicit, lambda, params)

%     ___________________________________________________
%
%                  THE WITF MODEL
%     ___________________________________________________
%


%% Set parameters
IsInit = false;
ShowFit = 10;
NumRep   = 5;
NumItInRep = 5;
InitMethod = 2;
MaxIt = 50;
ConvTol = 1e-7;
AbsErr = 1e-6;
MaxALSIT = 2;
MaxProIT = 5;
Regulizer = 2;
ProTol = [];
rnParam = [];
ShowOut = false;

I = size(X{1},1);
K = max(size(X));


if exist('params','var')
    if isfield(params,'initFac')
        IsInit = true;
        U = params.initFac{1};
        V = params.initFac{2};
        C = params.initFac{3};
        P = params.initFac{4};
    end
    
    % Show informantion
    if isfield(params,'ShowOut'), ShowOut = params.ShowOut; end
    
    % Show fit every 'ShowFit' iteration
    if isfield(params,'ShowFit'), ShowFit = params.ShowFit; end
    
    %Number of repetead initial analyses
    if isfield(params, 'NumRep'), NumRep = params.NumRep; end
    
    % Number of iterations in each initial fit
    if isfield(params, 'NumItInRep'), NumItInRep = params.NumItInRep; end
    
    % Max iteration
    if isfield(params, 'MaxIt'), MaxIt = params.MaxIt; end
    
    % Convergence Tol
    if isfield(params, 'Tol'), ConvTol = params.Tol; end
    
    % Fit error to exit
    if isfield(params, 'AbsErr'), AbsErr = params.AbsErr; end
    
    % Init Method: 0 = SVD+Random, 1 = SVD, else Random
    if isfield(params, 'InitMethod'), InitMethod = params.InitMethod; end
    
    % Max iteration of ALS Fitting
    if isfield(params, 'MaxALSIT'), MaxALSIT = params.MaxALSIT; end
    
    % Max iteration of Orthornormal Procrustes Fitting
    if isfield(params, 'MaxProIT'), MaxProIT = params.MaxProIT; end
    
    % Tol for WLSProcruster
    if isfield(params, 'ProTol'), ProTol = params.ProTol; end
    
    % Regulizer: 1 = L1 norm, else L2 norm
    if isfield(params, 'Regulizer'), Regulizer = params.Regulizer; end
    
    % number of noise for regularization
    if isfield(params, 'RandNoise')
        rnParam = parsenoiseparam(params.RandNoise);
    end
end

useRndNoise = ~isempty(rnParam) && rnParam.Tau > 0 && rnParam.Size > 0;

if nargin < 6 || isempty(lambda)
    lambdaU = 0;
    lambdaV = 0;
    lambdaC = 0;
else
    lambdaU = lambda(1);
    lambdaV = lambda(2);
    lambdaC = lambda(3);
end

if ~(length(size(X))==3||iscell(X))
    error('X must be a three-way array or a cell array')
end

LogInfo('\n\n Convergence criterion        : %g',ConvTol);
LogInfo('\n Maximal number of iterations : %d',MaxIt);
LogInfo('\n Number of factors            : %d\n',R);


% Initialize by ten small runs
if ~IsInit
    if InitMethod == 0
        LogInfo('\n Use best of %d initially fitted models', NumRep);
        
        paramInit = params;
        paramInit.MaxIt = NumItInRep;
        paramInit.InitMethod = 1;
        [U,V,C,P,bestfit]=WITF(X,R,W,wX,implicit,lambda,paramInit);
        for i = 2:NumRep
            paramInit.InitMethod = 2;
            
            LogInfo('\n Random initialization: %d', i);
            [Ur,Vr,Cr,Pr,fit]=WITF(X,R,W,wX,implicit,Opt,lambda,paramInit);
            if fit<bestfit
                U=Ur;V=Vr;C=Cr;P=Pr;
                bestfit = fit;
            end
        end
        % Initialize by SVD
    elseif InitMethod == 1
        LogInfo('\n SVD based initialization');
        
        XtX = zeros(I);
        sumwX = sum(wX);
        parfor k = 1:K
            meanXk = meandata(X{k}, implicit);
            Xk = full(X{k});
            Xk = Xk + meanXk - meanXk*logical(X{k});
            XtX = XtX + (wX(k)/sumwX)*(Xk*Xk');
        end
        [U,~,~] = svd(XtX,0);
        U = U(:,1:R);
        C = ones(K,R) + 0.1*randn(K,R);
        V = eye(R);
    else
        LogInfo('\n Random initialization')
        
        U = rand(I,R);
        C = rand(K,R);
        V = eye(R);
    end
    
    % init P and impute
    parfor k = 1:K
        meanXk = meandata(X{k}, implicit);
        Xk = full(X{k});
        Xk = Xk + meanXk - meanXk*logical(X{k});
        Qk = Xk' * (bsxfun(@times, U, C(k,:))*V');
        [Ur,~,Vr]  = svd(Qk,'econ');
        P{k}     = Ur*Vr';
    end
end

if useRndNoise
    LogInfo('\n\n Generating random noise\n');
    Wi = cell(1,K);
    wt = rnParam.Weight;
    if isscalar(wt), wt = wt.*ones(1,K); end
    for k = 1:K
        [X{k}, Wi{k}] = imputenoise(X{k}, W{k}, implicit, rnParam);
        W{k} = W{k} + wt(k)*Wi{k};
    end
end

wMax = ones(I,K);
parfor k = 1:K
     rowMax = max(W{k},[],2);
     rowMax(rowMax == 0) = eps;
     wMax(:,k) = rowMax;
end;

if isempty(ProTol), ProTol = zeros(1,K); end

HH = cell(K,I);
MSK = cell(K,I);
for k = 1:K
    Wk = wX(k)*W{k};
    parfor i=1:I
        HH{k,i} = Wk(i,:);
        MSK{k,i} = logical(HH{k,i});
    end
end

LogInfo('\n Fitting model ...')

fit = loss(X,W,wX,U,V,C,P,lambdaU,lambdaV,lambdaC,implicit,Regulizer);
Y = zeros(I,R,K);
tol = inf;
tid = tic;
it     = 0;
% Iterative part
while true
    oldfit = fit;
    it   = it + 1;
    
    parfor k = 1:K
        Y(:,:,k) = (wX(k).*W{k}.*X{k})*P{k};
    end
    
    % Update U,V,C
    [U,V,C]=alsfitting(Y,MSK,HH,wX,U,V,C,P,lambdaU,lambdaV,lambdaC,implicit,Regulizer,MaxALSIT);
    %[U,V,C] = cgfitting(Y,HH,U,V,C,P,lambdaU,lambdaV,lambdaC);
    
    % Update P
    P = WLSProcrustes( X, W, wMax, U, V, C, P, implicit, MaxProIT, ProTol);
    
    [fit] = loss(X,W,wX,U,V,C,P,lambdaU,lambdaV,lambdaC,implicit,Regulizer);
    
    tol = abs((fit-oldfit)/oldfit);
    % Print interim result
    if rem(it,ShowFit)==0||it == 1
        LogInfo('\n FitErr: %g\tIter: %g\tTol: %g\tExpire: %g secs',fit,it,tol,toc(tid));
    end
    
    if tol < ConvTol || it >= MaxIt || fit < AbsErr, break; end;
    
    if useRndNoise && rem(it,rnParam.ResplIt)==0
        parfor k = 1:K
            X{k}(Wi{k}) = 0;
            W{k}(Wi{k}) = 0;
            [X{k}, Wi{k}] = imputenoise(X{k}, W{k}, implicit, rnParam);
            W{k} = W{k} + wt(k)*Wi{k};
        end
    end
end

if rem(it,ShowFit)~=0 %Show final fit if not just shown
    LogInfo('%g        %g        %g        %g secs\n',fit,it,tol,toc(tid));
end

    function LogInfo(format,varargin)
        if ShowOut
            fprintf(format,varargin{:});
        end
    end

end

function P = WLSProcrustes( X, W, wMax, U, V, C, P, implicit, maxit, tol)
K = size(C,1);
if nargin < 9, maxit = 5; end
if nargin < 10, tol = -ones(1,K); end

parfor k=1:K
    wm = wMax(:,k);
    invwm = 1./wm;
    Wk = bsxfun(@times, W{k}, invwm);
    Z = bsxfun(@times, U, C(k,:))*V';
    tolk = tol(k);
    ls_p = inf;
    for i=1:maxit
        ZPk = Z*P{k}';
        D = ZPk - X{k};
        Y = ZPk - Wk.*D;
        if implicit
            Y = Y - bsxfun(@times, D, invwm);
        end
        Qk = Y' * bsxfun(@times, Z, wm);
        [u,~,v] = svd(Qk, 'econ');
        P{k} = u*v';
        
        if i < maxit && tolk > 0
            CPV = bsxfun(@times, P{k}*V, C(k,:))';
            ls = loss_k( X{k}, W{k}, U, CPV, implicit);
            if (ls_p - ls) < tolk, break; end
            ls_p = ls;
        end
    end
end
end

function fit = loss_k( Xk, Wk, U, CPV, implicit)
fit = 0;
parfor i = 1:size(U,1)
    Xki = Xk(i,:);
    Wki = Wk(i,:);
    Ui = U(i,:)*CPV;
    Xki2 = Xki.*Xki;
    Ui2 = Ui.*Ui;
    XkiUi = Xki.*Ui;
    lsi = sum(Xki2 .* Wki) + sum(Ui2 .* Wki) - 2*sum(XkiUi .* Wki);
    if implicit
        lsi = lsi + (sum(Xki2) + sum(Ui2) - 2*sum(XkiUi));
    end
    fit = fit + lsi;
end
end

function fit=loss(X,W,wX,U,V,C,P,lambdaU,lambdaV,lambdaC,implicit,Regulizer)
% Calculate fit and impute missing elements from model
fit = 0;
for k = 1:size(C,1)
    Wk = W{k};
    Xk = X{k};
    CPV = bsxfun(@times, P{k}*V, C(k,:))';
    fit = fit + wX(k)*loss_k( Xk, Wk, U, CPV, implicit);
end

if lambdaU~=0
    if Regulizer == 1
        fit = fit + lambdaU*sum(sum(abs(U)));
    else
        fit = fit + lambdaU*sum(sum(U.^2));
    end
end
if lambdaV~=0
    if Regulizer == 1
        fit = fit + lambdaV*sum(sum(abs(V)));
    else
        fit = fit + lambdaV*sum(sum(V.^2));
    end
end
if lambdaC~=0
    if Regulizer == 1
        fit = fit + lambdaC*sum(sum(abs(C)));
    else
        fit = fit + lambdaC*sum(sum(C.^2));
    end
end
end

function [U,V,C] = alsfitting(Y,MSK,HH,wX,U,V,C,P,lambdaU,lambdaV,lambdaC,implicit,Regulizer,maxit)
% Initialization
if nargin < 12, Regulizer = 2; end
if nargin < 13, maxit = 5; end

[I,J,K] = size(Y);
R = size(C,2);
regU = lambdaU*speye(R);
regV = lambdaV*speye(R);
regC = lambdaC*speye(R);

CTC = C'*C;
VTV = V'*V;
PVC = cell(1, K);
for it = 1:maxit
    
    % Update Ui
    
    CV = CTC.*VTV;
    XO = 0;
    
    % Khatri-Rao product
    PVCTPVC = 0;
    parfor k=1:K
        VC = bsxfun(@times,V,C(k,:));
        XO = XO + Y(:,:,k)*VC;
        PVC{k} = P{k}*VC;
        H = PVC{k}'*PVC{k};
        if implicit
            H = (1-wX(k))*H;
        end
        PVCTPVC = PVCTPVC + H;
    end
    
    parfor i=1:I
        H = CV;
        if lambdaU ~= 0
            if Regulizer == 1
                H = H + lambdaU*spdiags(sign(U(i,:))',0,R,R);
            else
                H = H + regU;
            end
        end
        HHi = HH(:,i);
        MSKi = MSK(:,i);
        for k=1:K
            MSKki = MSKi{k};
            HHki = HHi{k}(MSKki);
            PVCki = PVC{k}(MSKki,:);
            H = H +  bsxfun(@times,PVCki',HHki)*PVCki;
        end
        H = H - PVCTPVC;
        U(i,:) = XO(i,:)/H;
    end
    UTU = U'*U;
    
    % Update V
    H = CTC.*UTU;
    if lambdaV ~= 0
        if Regulizer == 1
            H = H + lambdaV*spdiags(sign(V(:)),0,R,R);
        else
            H = H + regV;
        end
    end
    H = kron(H,speye(R));
    
    XO = 0;
    % Khatri-Rao product
    parfor k=1:K
        XO = XO + Y(:,:,k)'*bsxfun(@times,U,C(k,:));
    end
    
    for k=1:K
        Ck = C(k,:);
        Pk = P{k};
        UCTUC = 0;
        MSKk = MSK(k,:);
        HHk = HH(k,:);
        parfor i=1:I
            UiCk = U(i,:).*Ck;
            UiCkTUiCk = UiCk'*UiCk;
            UCTUC = UCTUC + UiCkTUiCk;
            Pki = Pk(MSKk{i},:);
            HHki = HHk{i}(MSKk{i});
            H = H + kron(UiCkTUiCk, bsxfun(@times,Pki',HHki) * Pki);
        end
        PkTPk = Pk'*Pk;
        if implicit
            PkTPk = (1-wX(k))*PkTPk;
        end
        H = H - kron(UCTUC, PkTPk);
    end
    V = reshape(H\XO(:),R,R);
    VTV = V'*V;
    
    % Update Ck
    VU = VTV.*UTU;
    XO = 0;
    % Khatri-Rao product
    parfor j=1:J
        Xj = permute(Y(:,j,:),[3,1,2]);
        XO = XO + Xj * bsxfun(@times,U,V(j,:));
    end
    
    for k=1:K
        H = VU;
        if lambdaC ~= 0
            if Regulizer == 1
                regC = lambdaC*spdiags(sign(C(k,:))',0,R,R);
            end
            H = H + regC;
        end
        PV = P{k}*V;
        PVTPV = PV'*PV;
        if implicit
            PVTPV = (1-wX(k))*PVTPV;
        end
        HHk = HH(k,:);
        MSKk = MSK(k,:);
        parfor i=1:I
            Ui = U(i,:);
            PVU = bsxfun(@times,PV,Ui);
            PVUki = PVU(MSKk{i},:);
            HHki = HHk{i}(MSKk{i});
            H = H + ( bsxfun(@times,PVUki',HHki) * PVUki ...
                - bsxfun(@times, bsxfun(@times,PVTPV,Ui'), Ui) );
        end
        
        C(k,:) = XO(k,:)/H;
    end
    CTC = C'*C;
end
end
