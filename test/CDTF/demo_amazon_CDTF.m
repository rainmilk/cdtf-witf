load demo_amazon.mat

%% W-PARAFAC2
R = 25;

SN =Book_Matrix;

testIdx=missingInds_book_75;
missvalues = SN(testIdx);
SN(testIdx)=0;
X{1} = SN;

X{2} = Music_Matrix;

X{3} = DVD_Matrix;

X{4} = Video_Matrix;

Constraint = [0,0];
Opt = [0 0 2 0 0];

% Maximal number of iterations
Opt(2) = 10;

lambda = [0.05,0.001,0.001];
wX = [1,0.001,0.002,0.009];

epoch = 50;
MAE = inf;
[U,V,C,P,fit]=wparafac2(X,R,Constraint,Opt,lambda,wX);

for i=1:epoch
    params.initFac = {U,V,C,P};
    [U,V,C,P,fit]=wparafac2(X,R,Constraint,Opt,lambda,wX,params);

    predvalues = U*diag(sparse(C(1,:)))*(P{1}*V)';
    predvalues = predvalues(testIdx);
    predvalues( predvalues > 5 ) = 5;
    predvalues( predvalues < 1 ) = 1;
    MAE_i = sum(abs( missvalues - predvalues))/numel(missvalues);
    if MAE_i < MAE
        MAE = MAE_i;
    end
    fprintf('MAE: %g\n', MAE_i);
end

fprintf('Best MAE: %g\n', MAE);