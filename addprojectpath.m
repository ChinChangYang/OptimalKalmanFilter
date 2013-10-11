function addprojectpath
%ADDPROJECTPATH Add project pathes

% Get path information
currentPath = cd;
performancePath = sprintf('%s/performance', currentPath);
optimizerPath = sprintf('%s/robustoptimizer', currentPath);
solverPath = sprintf('%s/runner', currentPath);
utilPath = sprintf('%s/util', currentPath);

% Set path
addpath(genpath(performancePath));
addpath(genpath(optimizerPath));
addpath(genpath(solverPath));
addpath(genpath(utilPath));
end

