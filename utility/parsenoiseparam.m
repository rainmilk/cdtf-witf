function [ noiseParam ] = parsenoiseparam( noiseParam )
if ~isfield(noiseParam, 'ResplIt'), noiseParam.ResplIt = 1; end % regenerate noise every ResplIt iteration
if ~isfield(noiseParam, 'Tau'), noiseParam.Tau = inf; end % impute noisy data if # items < Tau
if ~isfield(noiseParam, 'Size'), noiseParam.Size = 10; end % the number of noisy data per user
if ~isfield(noiseParam, 'Std'), noiseParam.Std = .1; end % std of Gaussian noise
if ~isfield(noiseParam, 'Mean'), noiseParam.Mean = []; end % mean of Gaussian noise
if ~isfield(noiseParam, 'Weight'), noiseParam.Weight = .2*ones(1,K); end % weight on noise
end

