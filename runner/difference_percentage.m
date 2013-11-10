function difference_percentage
%DIFFERENCE_PERCENTAGE Maximal possible difference percentage between the
%perturbed function value and the function value in nominal case

load('minmin_results_by_jadebin_20131027.mat');
f0 = fmin;
load('maxmax_results_by_jadebin_20131027.mat');
f1 = -fmin;
p = abs(f0 - f1)/abs(f0);
fprintf('|maxmax - minmin|/|minmin| = %.2f%%\n', p * 100);
load('minmax_results_by_minmaxtcjadebin_201311081215.mat');
f1 = fminmax;
p = abs(f0 - f1)/abs(f0);
fprintf('|minmax - minmin|/|minmin| = %.2f%%\n', p * 100);

load('minminmin_results_by_jadebin_20131028.mat');
f0 = fmin;
load('maxmaxmax_results_by_jadebin_20131028.mat');
f1 = -fmin;
p = abs(f0 - f1)/abs(f0);
fprintf('|maxmaxmax - minminmin|/|minminmin| = %.2f%%\n', p * 100);
load('minminmax_results_by_minmaxtcjadebin_20131025.mat');
f1 = fminmax;
p = abs(f0 - f1)/abs(f0);
fprintf('|minminmax - minminmin|/|minminmin| = %.2f%%\n', p * 100);
load('maxminmax_results_by_maxminmaxtcjade_20131025.mat');
f1 = fbest;
p = abs(f0 - f1)/abs(f0);
fprintf('|maxminmax - minminmin|/|minminmin| = %.2f%%\n', p * 100);
end
