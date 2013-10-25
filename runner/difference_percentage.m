function [p1, p2] = difference_percentage
%DIFFERENCE_PERCENTAGE Maximal possible difference percentage between the
%perturbed function value and the function value in nominal case

load('minmax_results_by_minmaxtcjadebin.mat');
f0 = fminmax;
load('minminmax_results_by_minmaxtcjadebin.mat');
f1 = fminmax;
p1 = abs(f0 - f1)/abs(f0);
load('maxminmax_results_by_maxminmaxtcjade.mat');
f2 = fbest;
p2 = abs(f0 - f2)/abs(f0);

fprintf('|minminmax - minmax|/|minmax| = %.4E\n', p1);
fprintf('|maxminmax - minmax|/|minmax| = %.4E\n', p2);
end
