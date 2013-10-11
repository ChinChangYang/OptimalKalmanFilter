function display_inner_info(innerState)
%DISPLAY_INNER_INFO Display information of inner states

n_innerState = numel(innerState);
[D, ~] = size(innerState{1}.X);
stdX = zeros(D, n_innerState);
fstd = zeros(1, n_innerState);

for i = 1 : n_innerState
	stdX(:, i) = std(innerState{i}.X, 0, 2);
	fstd(i) = std(innerState{i}.f);
end

minStdX = min(stdX, [], 2);
maxStdX = max(stdX, [], 2);
meanStdX = mean(stdX(:));
medianStdX = median(stdX(:));

[maxmaxStdX, maxStdXIdx] = max(maxStdX);
[minminStdX, minStdXIdx] = min(minStdX);
fprintf(['InnerState: ', ...
	'Fstd: %0.2E,\t', ...
	'Min Xstd(%d): %0.2E,\t', ...
	'Max Xstd(%d): %0.2E,\t', ...
	'Mean Xstd: %0.2E,\t', ...
	'Median Xstd: %0.2E\n', ...
	'--\n'], ...
	mean(fstd), ...
	minStdXIdx, minminStdX, ...
	maxStdXIdx, maxmaxStdX, ...
	meanStdX, medianStdX);

end

