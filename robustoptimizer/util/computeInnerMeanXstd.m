function ret = computeInnerMeanXstd( innerState )
%COMPUTEINNERMEANXSTD Compute mean of standard deviation of solutions X in
%inner states

n_innerState = numel(innerState);
[D, ~] = size(innerState{1}.X);
stdX = zeros(D, n_innerState);

for i = 1 : n_innerState
	stdX(:, i) = std(innerState{i}.X, 0, 2);
end

ret = mean(stdX(:));
end

