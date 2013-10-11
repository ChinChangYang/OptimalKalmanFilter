function [xminmax1, xminmax2, fminmax, out] = minmaxdegl(fitfun, ...
	maxfunevals, lb1, ub1, lb2, ub2, options1, options2)
% MINMAXDEGL A differential evolution approach for solving constrained
% min-max optimization problems in Expert Systems with Applications 39
% (2013) 13440-13450 by Gilberto A.S. Segundo, Renato A. Krohling, and
% Rodrigo C. Cosme.
% MINMAXDEGL(fitfun, maxfunevals1, lb1, ub1, lb2, ub2) minimizes the
% function fitfun1 associated with a maximizer among box limitations [lb1,
% ub1] of minimizers and [lb2, ub2] of maximizers for the maximal function
% evaluations maxfunevals1.
% MINMAXDEGL(..., options1) minimizes the function with the given options
% Options1 for the 1st layer.
% MINMAXDEGL(..., options1, options2) minimize the function with the given
% options Options2 for the 2nd layer.

% Input arguments
if nargin <= 6
	options1 = [];
end

if nargin <= 7
	options2 = [];
end

D1 = numel(lb1);
D2 = numel(lb2);

% Default options for Layer 1
defaultOptions1.dimensionFactor = 5;
defaultOptions1.F = 0.7;
defaultOptions1.CR = 0.5;
defaultOptions1.maxStagnation = 20;
defaultOptions1.RecordPoint = 100;
defaultOptions1.TolX = eps;
defaultOptions1.TolFun = eps;
defaultOptions1.NeighborhoodRatio = 0.1;
options1 = setdefoptions(options1, defaultOptions1);
F = options1.F;
CR = options1.CR;
RecordPoint = options1.RecordPoint;
maxStagnation = options1.maxStagnation;
TolX = options1.TolX;
TolFun = options1.TolFun;
alpha = F;
beta = F;

% Default options for Layer 2
defaultOptions2 = defaultOptions1;
defaultOptions2.RecordPoint = 0;
options2 = setdefoptions(options2, defaultOptions2);

% Initialize algorithmic variables
NP1 = ceil(options1.dimensionFactor * D1);
nBest = round(0.2 * NP1);
isDisplayIter = strcmp(options1.Display, 'iter');

% Initialize variables
counteval = 0;
countiter = 1;
innerState = cell(1, NP1);
k = ceil(0.5 * (options1.NeighborhoodRatio * NP1));
w = 0.05 + 0.9 * rand(1, NP1);
wc = w;
out = initoutput(RecordPoint, D1, NP1, maxfunevals, ...
	'innerFstd', ...
	'innerMeanXstd');

% Initialize contour data
if isDisplayIter
	[XX, YY, ZZ] = minmaxcontourdata(D1, lb1, ub1, lb2, ub2, fitfun);
end

% Initialize population
X = zeros(D1, NP1);
for i = 1 : NP1
	X(:, i) = lb1 + (ub1 - lb1) .* rand(D1, 1);
end

% Evaluation
innerXbest = zeros(D2, NP1);
f = zeros(1, NP1);
for i = 1 : NP1
	[innerXbest(:, i), f(i), countevali, innerOut] = maxdegl(fitfun, lb2, ub2, X(:, i), [], options2);
	counteval = counteval + countevali;	
	innerState{i} = innerOut.final;
end

% Normalized fitness: Eq. (10)
nf = (f - min(f)) / (max(f) - min(f) + eps);

% Violations of lower and upper bounds
c_lb = zeros(D1, NP1);
c_ub = zeros(D1, NP1);
for i = 1 : NP1
	for j = 1 : D1
		c_lb(j, i) = max(0, lb1(j) - X(j, i));	% (13)
		c_ub(j, i) = max(0, X(j, i) - ub1(j));	% (13)
	end
end

% Cumulative sum of constraints violations: Eq. (17)
nvx = zeros(1, NP1);
vnew = zeros(1, NP1);
for i = 1 : NP1
	nvx(i) = sum(c_lb(:, i) > 0) + sum(c_ub(:, i) > 0);	
	
	for j = 1 : D1
		vnew(i) = vnew(i) + max(0, c_lb(j, i)) + max(0, c_ub(j, i));
	end
end

% Penalized fitness: Eq. (18)
vmax = max(vnew);
pf = zeros(1, NP1);
for i = 1 : NP1	
	if nvx(i) == 0
		pf(i) = nf(i);
	else
		pf(i) = nvx + vnew(i) / vmax;
	end
end

lastBestChange = countiter;

% Display
if isDisplayIter
	displayitermessages([X; innerXbest], [X; innerXbest], f, countiter, XX, YY, ZZ);
	display_inner_info(innerState);
end

% Record
out = updateoutput(out, X, f, counteval, ...
	'innerFstd', computeInnerFstd(innerState), ...
	'innerMeanXstd', computeInnerMeanXstd(innerState));

while true
	% Termination conditions
	outofmaxfunevals = counteval >= maxfunevals;
	fitnessconvergence = isConverged(f, TolFun);
	solutionconvergence = isConverged(X, TolX);
	
	% Convergence conditions	
	if outofmaxfunevals || fitnessconvergence || solutionconvergence
		break;
	end
	
	for i = 1 : NP1
		% Global best solution
		[~, g_best] = min(pf);
		
		% Neiborhoods index
		n_index = (i-k) : (i+k);
		lessthanone = n_index < 1;
		n_index(lessthanone) = n_index(lessthanone) + NP1;
		greaterthanNP = n_index > NP1;
		n_index(greaterthanNP) = n_index(greaterthanNP) - NP1;
		
		% Neiborhood solutions and fitness
		Xn = X(:, n_index);
		pfn = pf(n_index);
		
		% Best neiborhood
		[~, n_besti] = min(pfn);
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
		r1 = floor(1 + NP1 * rand);
		
		while i == r1
			r1 = floor(1 + NP1 * rand);
		end
		
		r2 = floor(1 + NP1 * rand);
		
		while i == r2 || r1 == r2
			r2 = floor(1 + NP1 * rand);
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
		
		xnew = wc(i) * gi + (1 - wc(i)) * Li;
		
		% Crossover
		jrand = floor(1 + rand * D1);
		for j = 1 : D1
			if rand >= CR && j ~= jrand
				xnew(j) = X(j, i);
			end
		end
		
		% Calculate the distance
		d = zeros(1, NP1);
		for j = 1 : NP1
			d(j) = norm(xnew - X(:, j));
		end
		
		% Sort distances
		[~, sort_d_index] = sort(d);
		bestPopC = innerXbest(:, sort_d_index(1 : nBest));
		
		% Optimize innerXbest
		[X2new, fnew, countevali, innerOutUi] = maxdegl(fitfun, lb2, ub2, xnew, bestPopC, options2);
		counteval = counteval + countevali;
		
		% Violations of lower and upper bounds
		c_lb_new = zeros(D1, 1);
		c_ub_new = zeros(D1, 1);
		for j = 1 : D1
			c_lb_new(j) = max(0, lb1(j) - xnew(j));	% (13)
			c_ub_new(j) = max(0, xnew(j) - ub1(j));	% (13)
		end
		
		% Cumulative sum of constraints violations: Eq. (17)
		nvx_new = sum(c_lb_new > 0) + sum(c_ub_new > 0);
		vnew_new = 0;
		
		for j = 1 : D1
			vnew_new = vnew_new + max(0, c_lb_new(j)) + max(0, c_ub_new(j));
		end
		
		% Selection
		if nvx_new == 0
			if fnew <= f(i)
				xnew_selected = true;
			else
				xnew_selected = false;
			end
		else
			if nvx_new < nvx(i)
				xnew_selected = true;
			elseif nvx_new == nvx(i) && vnew_new <= vnew(i)
				xnew_selected = true;
			else
				xnew_selected = false;
			end
		end
		
		if xnew_selected
			X(:, i) = xnew;
			innerXbest(:, i) = X2new;
			f(i) = fnew;
			nvx(i) = nvx_new;
			vnew(i) = vnew_new;
			c_lb(:, i) = c_lb_new;
			c_ub(:, i) = c_ub_new;	
			innerState{i} = innerOutUi.final;
			lastBestChange = countiter;

			% Normalized fitness: Eq. (10)
			nf = (f - min(f)) / (max(f) - min(f) + eps);
			
			% Penalized fitness: Eq. (18)
			vmax = max(vnew);
			for k = 1 : NP1
				if nvx(k) == 0
					pf(k) = nf(k);
				else
					pf(k) = nvx(k) + vnew(k) / vmax;
				end
			end
		end
	end
	
	% Display
	if isDisplayIter
		displayitermessages([X; innerXbest], [X; innerXbest], f, countiter, XX, YY, ZZ);		
		display_inner_info(innerState);
	end
	
	% Record
	out = updateoutput(out, X, f, counteval, ...
		'innerFstd', computeInnerFstd(innerState), ...
		'innerMeanXstd', computeInnerMeanXstd(innerState));
	
	countiter = countiter + 1;
end

% The best individual
[~, minpf_index] = min(pf);
xminmax1 = X(:, minpf_index);
xminmax2 = innerXbest(:, minpf_index);
fminmax = f(minpf_index);
out = finishoutput(out, X, f, counteval, ...
	'innerFstd', computeInnerFstd(innerState), ...
	'innerMeanXstd', computeInnerMeanXstd(innerState));
end

function [xmin, fmin, counteval, out] = maxdegl(fitfun, lb, ub, xnew, bestPopC, options)
% Input arguments
D = numel(lb);

% Default options
defaultOptions.dimensionFactor = 5;
defaultOptions.F = 0.7;
defaultOptions.CR = 0.5;
defaultOptions.maxStagnation = 20;
defaultOptions.maxIteration = 200;
defaultOptions.NeighborhoodRatio = 0.1;
options = setdefoptions(options, defaultOptions);
F = options.F;
CR = options.CR;
maxStagnation = options.maxStagnation;
maxIteration = options.maxIteration;
NP = ceil(options.dimensionFactor * D);
maxfunevals = maxIteration * NP;
alpha = F;
beta = F;

% Initialize variables
counteval = 0;
countiter = 1;
k = ceil(0.5 * (options.NeighborhoodRatio * NP));
w = 0.05 + 0.9 * rand(1, NP);
wc = w;
out = initoutput(0, D, NP, maxfunevals);

% Initialize Pop_B randomly
X = zeros(D, NP);
for i = 1 : NP
	X(:, i) = lb + (ub - lb) .* rand(D, 1);
end

% Copy bestPopC to Pop_B
if ~isempty(bestPopC)
	[~, nBest] = size(bestPopC);
	for i = 1 : nBest
		X(:, i) = bestPopC(:, i);
	end
end

% Evaluation
f = zeros(1, NP);
for i = 1 : NP
	f(i) = feval(fitfun, xnew, X(:, i));
	counteval = counteval + 1;
end

% Normalized fitness: Eq. (11)
nf = (max(f) - f) / (max(f) - min(f) + eps);

% Violations of lower and upper bounds
c_lb = zeros(D, NP);
c_ub = zeros(D, NP);
for i = 1 : NP
	for j = 1 : D
		c_lb(j, i) = max(0, lb(j) - X(j, i));	% (13)
		c_ub(j, i) = max(0, X(j, i) - ub(j));	% (13)
	end
end

% Cumulative sum of constraints violations: Eq. (17)
nvx = zeros(1, NP);
vnew = zeros(1, NP);
for i = 1 : NP
	nvx(i) = sum(c_lb(:, i) > 0) + sum(c_ub(:, i) > 0);	
	
	for j = 1 : D
		vnew(i) = vnew(i) + max(0, c_lb(j, i)) + max(0, c_ub(j, i));
	end
end

% Penalized fitness: Eq. (18)
vmax = max(vnew);
pf = zeros(1, NP);
for i = 1 : NP	
	if nvx(i) == 0
		pf(i) = nf(i);
	else
		pf(i) = nvx(i) + vnew(i) / vmax;
	end
end

lastBestChange = countiter;

% Record
out = updateoutput(out, X, f, counteval);

while countiter < maxIteration && countiter - lastBestChange < maxStagnation
	for i = 1 : NP
		% Global best solution
		[~, g_best] = min(pf);
		
		% Neiborhoods index
		n_index = (i-k) : (i+k);
		lessthanone = n_index < 1;
		n_index(lessthanone) = n_index(lessthanone) + NP;
		greaterthanNP = n_index > NP;
		n_index(greaterthanNP) = n_index(greaterthanNP) - NP;
		
		% Neiborhood solutions and fitness
		Xn = X(:, n_index);
		pfn = pf(n_index);
		
		% Best neiborhood
		[~, n_besti] = min(pfn);
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
		
		ynew = wc(i) * gi + (1 - wc(i)) * Li;
		
		% Crossover
		jrand = floor(1 + rand * D);
		for j = 1 : D
			if rand >= CR && j ~= jrand
				ynew(j) = X(j, i);
			end
		end
		
		% Evaluation
		fnew = feval(fitfun, xnew, ynew);
		counteval = counteval + 1;
		
		% Violations of lower and upper bounds
		c_lb_new = zeros(D, 1);
		c_ub_new = zeros(D, 1);
		for j = 1 : D
			c_lb_new(j) = max(0, lb(j) - ynew(j));	% (13)
			c_ub_new(j) = max(0, ynew(j) - ub(j));	% (13)
		end
		
		% Cumulative sum of constraints violations: Eq. (17)
		nvx_new = sum(c_lb_new > 0) + sum(c_ub_new > 0);
		vnew_new = 0;
		
		for j = 1 : D
			vnew_new = vnew_new + max(0, c_lb_new(j)) + max(0, c_ub_new(j));
		end
		
		% Selection		
		if nvx_new == 0
			if fnew >= f(i)
				ynew_selected = true;
			else
				ynew_selected = false;
			end
		else
			if nvx_new < nvx(i)
				ynew_selected = true;
			elseif nvx_new == nvx(i) && vnew_new <= vnew(i)
				ynew_selected = true;
			else
				ynew_selected = false;
			end
		end
		
		if ynew_selected
			X(:, i) = ynew;
			f(i) = fnew;
			nvx(i) = nvx_new;
			vnew(i) = vnew_new;
			c_lb(:, i) = c_lb_new;
			c_ub(:, i) = c_ub_new;	
			lastBestChange = countiter;		
			
			% Normalized fitness: Eq. (11)
			nf = (max(f) - f) / (max(f) - min(f) + eps);
			
			% Penalized fitness: Eq. (18)
			vmax = max(vnew);
			for k = 1 : NP
				if nvx(k) == 0
					pf(k) = nf(k);
				else
					pf(k) = nvx(k) + vnew(k) / vmax;
				end
			end
		end
	end
	
	countiter = countiter + 1;
	
	% Record
	out = updateoutput(out, X, f, counteval);
end

% The best individual
[~, minpf_index] = min(pf);
xmin = X(:, minpf_index);
fmin = f(minpf_index);

out = finishoutput(out, X, f, counteval);
end