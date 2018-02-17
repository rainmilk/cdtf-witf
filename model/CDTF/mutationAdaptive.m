function mutationChildren = mutationAdaptive(parents,options,GenomeLength,FitnessFcn,state,thisScore,thisPopulation,scale,shrink)
% MutationFcn
% Gaussian mutation

std = ( 1 - state.Generation/options.Generations);
mutationChildren = zeros(length(parents),GenomeLength);
for i=1:length(parents)
    parent = thisPopulation(parents(i),:);
    mutationChildren(i,:) = abs(parent  +  std * parent.*randn(1,length(parent)));
end
end
