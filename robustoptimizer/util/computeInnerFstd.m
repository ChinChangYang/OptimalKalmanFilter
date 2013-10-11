function ret = computeInnerFstd( innerState )
%COMPUTEINNERFSTD Compute standard deviation of function values in inner
%states

n_innerState = numel(innerState);
fstd = zeros(1, n_innerState);

for i = 1 : n_innerState
	fstd(i) = std(innerState{i}.f);
end

ret = mean(fstd);
end

