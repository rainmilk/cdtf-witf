p = mfilename('fullpath');
dir=fileparts(p);
base = [dir filesep];

addpath([base 'metrics']);
addpath([base 'data']);

addpath([base 'model']);
addpath([base 'model/CDTF']);
addpath([base 'model/MF']);

addpath([base 'test']);
addpath([base 'test/CDTF']);

addpath([base 'utility']);

clear




