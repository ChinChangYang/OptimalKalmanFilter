function verify_minmax
%VERIFY_MINMAX Verify the minmax solution

fminmax = 0.0120347390864446;
xminmax1 = [0.499994432810761;0.041377184011307;0.249768724165629;0.900005183763698];
xminmax2 = [0.500000000000000;-0.050000000000000;-0.050000000000000;1.100000000000000];
counter = 0;

% Verify f(xminmax1, xminmax2)
f = error_nozeta(xminmax1, xminmax2);
if abs(fminmax - f) > 1e4 * abs(f)
	counter = counter + 1;
	fprintf('Failed: f(xminmax1, xminmax2) ~= fminmax\n');
end

% Verify max(f(xminmax1, x))
fitfun = @(x) -error_nozeta(xminmax1, x);
maxfunevals = 1e6;
lb = [0.3; -0.05; -0.05; 0.9];
ub = [0.5; 0.25; 0.45; 1.1];
solverOptions.dimensionFactor = 5;
[xmin, fmin, ~] = ...
	feval('jadebin', fitfun, lb, ub, maxfunevals, solverOptions);

if norm(xmin - xminmax2) > 1e-4
	counter = counter + 1;
	fprintf('Failed: xminmax2 values is not consistent\n');
end

% Verify xminmax1
x_hat = xminmax1 .* (1 + 1e-4 * randn(4, 1));
for i = 1 : numel(lb)
	if x_hat(i) < lb(i)
		x_hat(i) = 2 * lb(i) - x_hat(i);
	elseif x_hat(i) > ub(i)
		x_hat(i) = 2 * ub(i) - x_hat(i);
	end
	
	if x_hat(i) < lb(i)
		x_hat(i) = lb(i);
	elseif x_hat(i) > ub(i)
		x_hat(i) = ub(i);
	end
end

fitfun = @(x) -error_nozeta(x_hat, x);
maxfunevals = 1e6;
lb = [0.3; -0.05; -0.05; 0.9];
ub = [0.5; 0.25; 0.45; 1.1];
[~, fmin, ~] = ...
	feval('jadebin', fitfun, lb, ub, maxfunevals, solverOptions);

if abs(fminmax + fmin) > 1e-7
	counter = counter + 1;	
	fprintf('Failed: xminmax1 is not the minimaximizer\n');
end

if counter == 0
	fprintf('All Pass!\n');
else
	fprintf('Failed: %d\n', counter);
end
end

