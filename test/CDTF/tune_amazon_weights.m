%% W-PARAFAC2
R = 25;
domain='book';
percent=0.75;
SN =Book_Matrix;%ratings(1:2000,:)';
testIdx=missingInds_book_75;
missvalues = SN(testIdx);
SN(testIdx)=0;
X{1}=SN;

testIdx = sparse(testIdx);

SN=DVD_Matrix;
X{3}=SN;


SN=Video_Matrix;
X{4}=SN;

Constraint = [0,0];
Opt = [0 0 0 0 0];
lambda = [0.05,0.001,0.001];

% Coordinate Ascend Learning to tune weights
K = length(X);
% wX = [1, 0.2*ones(1,K-1)];
% [A,H,C,P,fit]=wparafac2(X,R,Constraint,Opt,lambda,wX);
tol = 5e-4;
% Maximal number of iterations
Opt(2) = 30;
cache = containers.Map('KeyType','char','ValueType','double');

fGA =  @(weight) fMAE( weight, cache, R, lambda, X, testIdx, missvalues, tol, Constraint, Opt);
p1 = diag([0.01,0.1,1,10])*ones(4,K-1);
initPop = [1/3*p1;2/3*p1;p1];
options = gaoptimset('StallGenLimit',5,...
    'PopulationSize',length(initPop),'PlotFcns',@gaplotbestf,'InitialPopulation',initPop,...
    'MutationFcn',@mutationAdaptive);
[x,fval,exitflag,output,population,scores] = ga(fGA,K-1,[],[],[],[],[],[],[],options);

fprintf('Best weights %s\n', mat2str(x));
