function addprojectpath
%ADDPROJECTPATH Add project pathes

% Get path information
currentPath = cd;
optimizerPath = sprintf('%s/robustoptimizer', currentPath);
utilPath = sprintf('%s/util', currentPath);

% Problem path
problemPath = sprintf('%s/Problem1', currentPath);
% problemPath = sprintf('%s/Problem3', currentPath);

% Set path
addpath(genpath(problemPath));
addpath(genpath(optimizerPath));
addpath(genpath(utilPath));
end

