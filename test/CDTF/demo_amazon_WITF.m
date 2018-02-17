load demo_amazon.mat

%% W-PARAFAC2
R = 25;

TargetDomain = 1;

implicit = false;

testIdx = missingInds_book_75;

TD = Book_Matrix;
[i,j] = find(testIdx);
testMat = sparse(i,j,TD(testIdx),size(TD,1),size(TD,2));
TD(testIdx) = 0;


Y{1} = TD;

Y{2} = Music_Matrix;

Y{3} = DVD_Matrix;

Y{4} = Video_Matrix;

W = cell(size(Y));
for k = 1:length(W)
    W{k} = logical(Y{k});
end

%% WITF
alpha = .01;
nO = cellfun(@nnz, Y);
nAux = nO(2:end);
wX = nO(1)./ nAux;
wX = [1, alpha .* wX];

lambda = 1*ones(1,3);

params = struct;
params.ShowOut = 1;
params.InitMethod = 2;
params.MaxIt = 3;
params.MaxALSIT = 1;
params.MaxProIT = 1;
params.ShowFit = 1;

% params.RandNoise.Tau = 20;
% params.RandNoise.Weight = 0.2;
% params.RandNoise.Size = 10;
% params.RandNoise.Std = .05;
% params.RandNoise.ResplIt = 1;

[U,V,C,P,fit]=WITF(Y,R,W,wX,implicit,lambda,params);

params.MaxIt = 3;

epoch = 50;
for i=1:epoch
    %params.ImputeParam.Wi = Wi;  
    params.initFac = {U,V,C,P};
    [U,V,C,P,fit]=WITF(Y,R,W,wX,implicit,lambda,params);
    UT = bsxfun(@times, U, C(1,:));
    VT = P{1}*V;
    %MAE_TR = MAE_RMSE(UT', VT', TD);
    [MAE_witf,RMSE_witf] = MAE_RMSE(UT', VT', testMat);
    fprintf('\nWITF: Round: %d, R: %d, MAE: %g, RMSE: %g', i, R, MAE_witf, RMSE_witf);
end

%% WRMF
maxit = 5;
epoch = 20;

csuser = ~any(W{1},2);
Uold = UT(csuser,:);

muU = UT;
muV = VT;

U_FT = muU + 0.01*randn(size(UT));
V_FT = muV + 0.01*randn(size(VT));
lambda_U = 0.1;
lambda_V = 0.1;

for i=1:epoch  
    [ U_FT, V_FT ] = WRMF( Y{1}, W{1}, lambda_U, lambda_V,...
        U_FT, V_FT, muU, muV, 'MaxIt', maxit, 'CplxCtrl', 1);
    U_C = U_FT;
    U_C(csuser,:) = Uold;
    %U_C(csuser,:) = MF_Rotation(Uold, VT, V_FT);
    [MAE_wrmf,RMSE_wrmf] = MAE_RMSE(U_C', V_FT', testMat);
    fprintf('\nFine Tuning: Round: %d, R: %d, MAE: %g, RMSE: %g', i, R, MAE_wrmf, RMSE_wrmf);
end