%% Reference to algorithms
%   Liang Hu, Jian Cao, Guandong Xu, Longbing Cao, Zhiping Gu, Can Zhu
%   Personalized recommendation via cross-domain triadic factorization. WWW 2013: 595-606
%  Liang Hu, Longbing Cao, Jian Cao, Zhiping Gu, Guandong Xu, Dingyu Yang:
Learning Informative Priors from Heterogeneous Domains to Improve Recommendation in Cold-Start User Domains. ACM Trans. Inf. Syst. 35(2): 13:1-13:37 (2016)

%% initpath.m
%	Initialize the paths (Running it first before calling other script)


%% Directory: model
%	Code for models

%% Directory: model/CDTF
%   CDTF related algorithms
%
% File List: 
%   wparafac2.m
%       Core algorithm for CDTF model
%
%   fMAE.m
%       Fitness function for GA
%
%   fAUC.m
%       Fitness function for GA
%
%   mutationAdaptive.m
%       MutationFcn for GA
%
%   if_input.m
%       Generate modified input for CDTF-IF
%
%   WITF.m
%       A refined implement of CDTF, which is published in TOIS. A GPU-version is desired.



%% Directory: test
% Code for testing models

%% Directory: test/CDTF
%   CDTF related algorithms

% File List: 

%  demo_amazon_CDTF.m
%      A demo to run CDTF algorithm over Amazon dataset with 75% traning data
%
%  demo_amazon_WITF.m
%      A demo to run WITF algorithm over Amazon dataset with 75% traning data
%
%  tune_amazon_weights.m
%      Find optimal weights using GA for Amazon data

%% Directory: metrics
%	Code for various evaluation metrics

%% Directory: data
%   Data files