function ret = isConverged(X, TolX)
stdX = std(X, 0, 2);
meanX = mean(abs(X), 2);
ret = all(stdX <= meanX * 100 * eps) || ...
	all(stdX <= 100 * realmin) || all(stdX <= TolX);
end
