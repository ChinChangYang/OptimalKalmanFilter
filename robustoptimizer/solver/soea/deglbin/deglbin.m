function [xmin, fmin, out] = deglbin(fitfun, lb, ub, maxfunevals, options)
% DEGLBIN Differential Evolution with global and local neighborhoods
% DEGLBIN(fitfun, lb, ub, maxfunevals) minimize the function fitfun in
% box constraints [lb, ub] with the maximal function evaluations
% maxfunevals.
% DEGLBIN(..., options) minimize the function by solver options.
if nargin <= 4
	options = [];
end

defaultOptions.dimensionFactor = 10;
defaultOptions.CR = 0.5;
defaultOptions.F = 0.7;
defaultOptions.NeighborhoodRatio = 0.1;
defaultOptions.Display = 'off';
defaultOptions.RecordPoint = 100;
defaultOptions.ftarget = -Inf;
defaultOptions.TolFun = eps;
defaultOptions.TolX = 100 * eps;
defaultOptions.initial.X = [];
defaultOptions.initial.f = [];
defaultOptions.initial.w = [];

options = setdefoptions(options, defaultOptions);
dimensionFactor = options.dimensionFactor;
CR = options.CR;
F = options.F;
isDisplayIter = strcmp(options.Display, 'iter');
RecordPoint = max(1, floor(options.RecordPoint));
ftarget = options.ftarget;
TolFun = options.TolFun;
TolX = options.TolX;
X = options.initial.X;
f = options.initial.f;
w = options.initial.w;

alpha = F;
beta = F;
D = numel(lb);

if isempty(X)
	NP = ceil(dimensionFactor * D);
else
	[~, NP] = size(X);
end

% Initialize variables
counteval = 0;
countiter = 1;
out = initoutput(RecordPoint, D, NP, maxfunevals);

% Initialize contour data
if isDisplayIter
	[XX, YY, ZZ] = preparecontourdata(D, lb, ub, fitfun);
end

% Initialize population
if isempty(X)
	if NP < 1e1
		LHS = lhsdesign(NP, D, 'iteration', 10)';
	elseif NP < 1e2
		LHS = lhsdesign(NP, D, 'iteration', 2)';
	else
		LHS = rand(D, NP);
	end
	
	X = zeros(D, NP);
	for i = 1 : NP
		X(:, i) = lb + (ub - lb) .* LHS(:, i);
	end
end

% Evaluation
if isempty(f)
	f = zeros(1, NP);
	for i = 1 : NP
		f(i) = feval(fitfun, X(:, i));
		counteval = counteval + 1;
	end
end

% Initialize variables
if isempty(w)	
	w = 0.05 + 0.9 * rand(1, NP);
end

k = ceil(0.5 * (options.NeighborhoodRatio * NP));
wc = w;
V = X;
U = X;

% Display
if isDisplayIter
	displayitermessages(X, U, f, countiter, XX, YY, ZZ);
end

% Record
out = updateoutput(out, X, f, counteval);

% Iteration counter
countiter = countiter + 1;

while true
	% Termination conditions
	stdf = std(f);
	stdX = std(X, 0, 2);
	meanX = mean(X, 2);
	outofmaxfunevals = counteval > maxfunevals - NP;
	reachftarget = min(f) <= ftarget;
	fitnessconvergence = stdf <= mean(abs(f)) * 100 * eps || stdf <= realmin || stdf <= TolFun;
	solutionconvergence = all(stdX <= meanX * 100 * eps) || all(stdX <= 100 * realmin) || ...
		all(stdX <= TolX);
	
	% Convergence conditions
	if outofmaxfunevals || reachftarget || fitnessconvergence || solutionconvergence
		break;
	end
	
	% Mutation
	% Global best
	[~, g_best] = min(f);
	
	for i = 1 : NP
		% Neiborhoods index
		n_index = (i-k) : (i+k);
		lessthanone = n_index < 1;
		n_index(lessthanone) = n_index(lessthanone) + NP;
		greaterthanNP = n_index > NP;
		n_index(greaterthanNP) = n_index(greaterthanNP) - NP;
		
		% Neiborhood solutions and fitness
		Xn = X(:, n_index);
		fn = f(n_index);
		
		% Best neiborhood
		[~, n_besti] = min(fn);
		Xn_besti = Xn(:, n_besti);
		
		% Random neiborhood index
		n_index(n_index == i) = [];
		Xn = X(:, n_index);
		p = ceil(rand * numel(n_index));
		q = ceil(rand * numel(n_index));
		
		while p == q
			q = ceil(rand * numel(n_index));
		end
		
		% Random neiborhood solutions
		Xp = Xn(:, p);
		Xq = Xn(:, q);
		
		% Local donor vector
		Li = X(:, i) + alpha * (Xn_besti - X(:, i)) + ...
			beta * (Xp - Xq);
		
		% Global donor vector
		r1 = floor(1 + NP * rand);
		
		while i == r1
			r1 = floor(1 + NP * rand);
		end
		
		r2 = floor(1 + NP * rand);
		
		while i == r2 || r1 == r2
			r2 = floor(1 + NP * rand);
		end
		
		gi = X(:, i) + alpha * (X(:, g_best) - X(:, i)) + ...
			beta * (X(:, r1) - X(:, r2));
		
		% Self-adaptive weight factor
		wc(i) = w(i) + F * (w(g_best) - w(i)) + ...
			F * (w(r1) - w(r2));
		
		if wc(i) < 0.05
			wc(i) = 0.05;
		elseif wc(i) > 0.95
			wc(i) = 0.95;
		end
		
		V(:, i) = wc(i) * gi + (1 - wc(i)) * Li;
	end
	
	for i = 1 : NP
		% Binominal Crossover
		jrand = floor(1 + D * rand);
		for j = 1 : D
			if rand < CR || j == jrand
				U(j, i) = V(j, i);
			else
				U(j, i) = X(j, i);
			end
		end
	end
	
	% Repair
	for i = 1 : NP
		for j = 1 : D
			if U(j, i) < lb(j)
				U(j, i) = X(j, i) + rand * (lb(j) - X(j, i));
			elseif U(j, i) > ub(j)
				U(j, i) = X(j, i) + rand * (ub(j) - X(j, i));
			end
		end
	end
	
	% Display
	if isDisplayIter
		displayitermessages(X, U, f, countiter, XX, YY, ZZ);
	end
	
	% Selection
	counteval = counteval + NP;
	for i = 1 : NP
		fui = feval(fitfun, U(:, i));
		
		if fui < f(i)
			f(i) = fui;
			X(:, i) = U(:, i);
			w(i) = wc(i);
		end
	end
	
	% Record
	out = updateoutput(out, X, f, counteval);
	
	% Iteration counter
	countiter = countiter + 1;
end

[fmin, fminidx] = min(f);
xmin = X(:, fminidx);

if fmin < out.bestever.fmin
	out.bestever.fmin = fmin;
	out.bestever.xmin = xmin;
end

% Final state of this algorithm
final.w = w;
out = finishoutput(out, X, f, counteval, 'final', final);
end