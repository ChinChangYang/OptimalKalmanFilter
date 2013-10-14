function [xbest1, xbest2, xbest3, fbest, out] = maxminmaxtcjade(...
	fitfun, ...		% objective function f(x,y,z)
	maxfunevals, ...	% maximal function evaluations
	lb1, ub1, ...		% lower and upper bounds for the 1st layer
	lb2, ub2, ...		% lower and upper bounds for the 2nd layer
	lb3, ub3, ...		% lower and upper bounds for the 3rd layer
	options1, ...		% options for the 1st layer
	options2, ...		% options for the 2nd layer
	options3)			% options for the 3rd layer
% MAXMINMAXTCJADE Max-Min-Max Tracer JADE
% [xbest1, xbest2, xbest3] = MAXMINMAXTCJADE(fitfun, maxfunevals, lb1,
% ub1, lb2, ub2, lb3, ub3) maximizes the function fitfun associated with a
% maximizer xbest1 among box limitations [lb1, ub1], a minimizer xbest2
% among [lb2, ub2], and a maximizer xbest3 among [lb3, ub3], which are
% searched by evolutionary algorithm within a maximal function evaluations
% maxfunevals.
% MAXMINMAXTCJADE(..., options1) maximizes the function with the given
% options Options1 for the 1st layer.
% MAXMINMAXTCJADE(..., options1, options2) maximizes the function with
% the given options Options2 for the 2nd layer.
% MAXMINMAXTCJADE(..., options1, options2, options3) maximizes the
% function with the given options Options3 for the 3rd layer.
% [..., fbest] = MAXMINMAXTCJADE(...) returns the function value of the
% max-min-max solution.
if nargin <= 9
	options1 = [];
end

if nargin <= 10
	options2 = [];
end

if nargin <= 11
	options3 = [];
end

D1 = numel(lb1);
D2 = numel(lb2);
D3 = numel(lb3);

% Default options for Layer 1
defaultOptions1.dimensionFactor = 10;
defaultOptions1.F = 0.9;
defaultOptions1.CR = 0.99;
defaultOptions1.delta_F = 0.1;
defaultOptions1.delta_CR = 0.01;
defaultOptions1.p = 0.05;
defaultOptions1.w = 0.1;
defaultOptions1.Display = 'off';
defaultOptions1.RecordPoint = 100;
defaultOptions1.TolX = 0;
defaultOptions1.TolFun = 0;
defaultOptions1.TolStagnationIteration = 20;
defaultOptions1.TolX_DecayRate = 0.5;
defaultOptions1.InnerSolver = 'minmaxtcjadebin';
defaultOptions1.initial.X = [];
defaultOptions1.initial.f = [];
defaultOptions1.initial.mu_F = [];
defaultOptions1.initial.mu_CR = [];
defaultOptions1.initial.innerState = [];
options1 = setdefoptions(options1, defaultOptions1);

% Default options for Layer 2
defaultOptions2.dimensionFactor = 10;
defaultOptions2.Display = 'off';
defaultOptions2.RecordPoint = 0;
defaultOptions2.TolFun = 0;
defaultOptions2.TolX = 0;
defaultOptions2.InnerSolver = 'jadebinwoa';
options2 = setdefoptions(options2, defaultOptions2);

% Default options for Layer 3
defaultOptions3.dimensionFactor = 10;
options3 = setdefoptions(options3, defaultOptions3);

% Initialize algorithmic variables
dimensionFactor = max(1, options1.dimensionFactor);
delta_F = options1.delta_F;
delta_CR = options1.delta_CR;
p = options1.p;
w = options1.w;
isDisplayIter = strcmp(options1.Display, 'iter');
RecordPoint = max(0, floor(options1.RecordPoint));
TolFun = options1.TolFun;
TolX = options1.TolX;
TolStagnationIteration = options1.TolStagnationIteration;
TolX_DecayRate = options1.TolX_DecayRate;
innerSolver = options1.InnerSolver;
X = options1.initial.X;
f = options1.initial.f;
mu_F = options1.initial.mu_F;
mu_CR = options1.initial.mu_CR;
innerState = options1.initial.innerState;
existInnerState = ~isempty(innerState);

NP1 = ceil(dimensionFactor * D1);
NP2 = ceil(options2.dimensionFactor * D2);
NP3 = ceil(options3.dimensionFactor * D3);

% Initialize contour data
if isDisplayIter
	plotFitfun = @(x, y) feval(fitfun, x, y, 0.5 * (lb3 + ub3));
	[XX, YY, ZZ] = minmaxcontourdata(D1, lb1, ub1, lb2, ub2, plotFitfun);
end

% Initialize population
if isempty(X)
	if NP1 < 1e1
		LHS = lhsdesign(NP1, D1, 'iteration', 10)';
	elseif NP1 < 1e2
		LHS = lhsdesign(NP1, D1, 'iteration', 2)';
	else
		LHS = rand(D1, NP1);
	end
	
	X = zeros(D1, NP1);
	for i = 1 : NP1
		X(:, i) = lb1 + (ub1 - lb1) .* LHS(:, i);
	end
end

% Initialize inner states
if isempty(innerState)
	innerState = cell(1, NP1);
end

% Initialize variables
counteval = 0;
countiter = 1;
countStagnation = 0;
successRate = 0;
X_Converged_FEs = zeros(1, NP1);
U_Converged_FEs = zeros(1, NP1);
innerState_Xstd = zeros(1, NP1);
innerXbest1 = zeros(D2, NP1);
innerXbest2 = zeros(D3, NP1);
innerUbest1 = innerXbest1;
innerUbest2 = innerXbest2;
V = X;
U = X;
V2 = zeros(D2, NP2, NP1);
V3 = zeros(D3, NP3, NP2, NP1);
pbest_size = p * NP1;
out = initoutput(RecordPoint, D1, NP1, maxfunevals, ...
	'innerFstd', ...
	'innerMeanXstd', ...
	'successRate', ...
	'X_Converged_FEs', ...
	'U_Converged_FEs', ...
	'mu_F', ...
	'mu_CR');

% Evaluation
if isempty(f)
	f = zeros(1, NP1);	
	innerMaxfunevalsX = NP2 * NP3;
	
	for i = 1 : NP1
		innerFitfun = @(y, z) feval(fitfun, X(:, i), y, z);	
		optionsX2i = options2;
		
		if existInnerState
			optionsX2i.initial = innerState{i};
			optionsX2i.initial.A = []; % disabled external option archive
		end
		
		[innerXbest1(:, i), innerXbest2(:, i), innerFbest, innerOut] = ...
			feval(innerSolver, innerFitfun, innerMaxfunevalsX, lb2, ub2, ...
			lb3, ub3, optionsX2i, options3);
		
		counteval = counteval + innerOut.fes(end);
		f(i) = -innerFbest;
		innerState{i} = innerOut.final;
	end
end

% Sort
[f, fidx] = sort(f);
X = X(:, fidx);
innerState = innerState(fidx);

% mu_F
if isempty(mu_F)
	mu_F = options1.F;
end

% mu_CR
if isempty(mu_CR)
	mu_CR = options1.CR;
end

% Display
if isDisplayIter
	displayitermessages([X; innerXbest1], [X; innerXbest1], f, countiter, ...
		XX, YY, ZZ, 'counteval', counteval, ...
		'successRate', successRate, ...
		'mu_F', mu_F, ...
		'mu_CR', mu_CR);
	
	display_inner_info(innerState);
end

% Record minimal function values
out = updateoutput(out, X, f, counteval, ...
	'innerFstd', computeInnerFstd(innerState), ...
	'innerMeanXstd', computeInnerMeanXstd(innerState), ...
	'successRate', successRate, ...
	'X_Converged_FEs', mean(X_Converged_FEs), ...
	'U_Converged_FEs', mean(U_Converged_FEs), ...
	'mu_F', mu_F, ...
	'mu_CR', mu_CR);

countiter = countiter + 1;

while true
	% Termination conditions
	outofmaxfunevals = counteval >= maxfunevals;
	fitnessconvergence = isConverged(f, TolFun);
	solutionconvergence = isConverged(X, TolX);
	stagnation = countStagnation >= TolStagnationIteration;
	
	% Convergence conditions
	if outofmaxfunevals || fitnessconvergence || solutionconvergence ...
			|| stagnation
		break;
	end
	
	% Scaling factor and crossover rate
	S_CR = zeros(1, NP1);
	CR = mu_CR + delta_CR * randn(1, NP1);
	CR(CR > 1) = 1 - eps;
	CR(CR < 0) = eps;
	S_F = zeros(1, NP1);
	F = cauchyrnd(mu_F, delta_F, NP1, 1);
	F(F > 1) = 1 - eps - 0.01 * rand;
	
	for retry = 1 : 3
		if all(F > 0)
			break;
		end
		
		F(F <= 0) = cauchyrnd(mu_F, delta_F, sum(F <= 0), 1);
		F(F > 1) = 1 - eps - 0.01 * rand;
	end
	
	F(F <= 0) = 0.01 * mu_F * (1 + rand);
	
	Succ_Counter = 0;
	
	% Mutation
	for i = 1 : NP1
		% Try generating V within bounds
		for retry_within_bounds = 1 : NP1
			
			% Generate pbest_idx
			for retry = 1 : 3
				pbest_idx = max(1, ceil(rand * pbest_size));
				if ~all(X(:, pbest_idx) == X(:, i))
					break;
				end
			end
			
			% Generate r1
			for retry = 1 : NP1
				r1 = floor(1 + NP1 * rand);
				if i ~= r1
					break;
				end
			end
			
			% Generate r2
			for retry = 1 : NP1 * NP1
				r2 = floor(1 + NP1 * rand);
				if ~(all(X(:, i) == X(:, r2)) || all(X(:, r1) == X(:, r2)))
					break;
				end
			end
			
			% Generate Vi
			V(:, i) = X(:, i) + F(i) .* ...
				(X(:, pbest_idx) - X(:, i) + X(:, r1) - X(:, r2));
			
			% Check boundary			
			if all(V(:, i) >= lb1) && all(V(:, i) <= ub1)
				
				% Prediction
				V2(:, :, i) = innerState{i}.X + F(i) .* ...
					(innerState{pbest_idx}.X - innerState{i}.X + ...
					innerState{r1}.X - innerState{r2}.X);
				
				for j = 1 : NP2
					V3(:, :, j, i) = innerState{i}.innerState{j}.X + ...
						F(i) .* ...
					(innerState{pbest_idx}.innerState{j}.X - ...
					innerState{i}.innerState{j}.X + ...
					innerState{r1}.innerState{j}.X - ...
					innerState{r2}.innerState{j}.X);
				end
				
				% Perturbation
				V2(:, :, i) = V2(:, :, i) .* (1 + 100 * eps * randn(D2, NP2));
				
				V3(:, :, :, i) = ...
					V3(:, :, :, i) .* (1 + 100 * eps * randn(D3, NP3, NP2));
				
				break;
			end
		end
	end
    
	% Crossover
	for i = 1 : NP1
		jrand = floor(1 + D1 * rand);
		
		for j = 1 : D1
			if rand < CR(i) || j == jrand
				U(j, i) = V(j, i);
			else
				U(j, i) = X(j, i);
			end
		end
	end
	
	% Repair
	for i = 1 : NP1
		for j = 1 : D1
			if U(j, i) < lb1(j)
				U(j, i) = 2 * lb1(j) - U(j, i);
			elseif U(j, i) > ub1(j)
				U(j, i) = 2 * ub1(j) - U(j, i);
			else
				continue;
			end
			
			if U(j, i) < lb1(j)
				U(j, i) = lb1(j);
			elseif U(j, i) > ub1(j)
				U(j, i) = ub1(j);
			end
		end
		
		for j = 1 : NP2
			for k = 1 : D2
				if V2(k, j, i) < lb2(k)
					V2(k, j, i) = 2 * lb2(k) - V2(k, j, i);
				elseif V2(k, j, i) > ub2(k)
					V2(k, j, i) = 2 * ub2(k) - V2(k, j, i);
				else
					continue;
				end
				
				if V2(k, j, i) < lb2(k)
					V2(k, j, i) = lb2(k);
				elseif V2(k, j, i) > ub2(k)
					V2(k, j, i) = ub2(k);
				end
			end
		end
		
		for j = 1 : NP2
			for k = 1 : NP3
				for m = 1 : D3
					if V3(m, k, j, i) < lb3(m)
						V3(m, k, j, i) = 2 * lb3(m) - V3(m, k, j, i);
					elseif V3(m, k, j, i) > ub3(m)
						V3(m, k, j, i) = 2 * ub3(m) - V3(m, k, j, i);
					else
						continue;
					end
					
					if V3(m, k, j, i) < lb3(m)
						V3(m, k, j, i) = lb3(m);
					elseif V3(m, k, j, i) > ub3(m)
						V3(m, k, j, i) = ub3(m);
					end					
				end
			end
		end
	end
	
	% TolX Decision
	for i = 1 : NP1
		innerState_Xstd(i) = mean(std(innerState{i}.X, [], 2));
	end
	
	innerTolX = TolX_DecayRate * mean(innerState_Xstd);
	
	% Selection
	successRate = 0;
	FailedIteration = true;
	innerMaxfunevalsX = 400 * NP2 * NP3;
	for i = 1 : NP1
		% Compute fxi, f(i)
		innerFitfunXi = @(y, z) feval(fitfun, X(:, i), y, z);
		optionsX2i = options2;
		optionsX2i.initial = innerState{i};
		optionsX2i.initial.A = []; % disabled external option archive
		optionsX2i.TolX = max(options2.TolX, innerTolX);
		
		[innerXbest1(:, i), innerXbest2(:, i), innerFbest, innerOutXi] = ...
			feval(innerSolver, innerFitfunXi, innerMaxfunevalsX, ...
			lb2, ub2, ...
			lb3, ub3, ...
			optionsX2i, options3);
		
		X_Converged_FEs(i) = innerOutXi.fes(end);
		counteval = counteval + innerOutXi.fes(end);
		f(i) = -innerFbest;
		
		% Compute fui
		innerFitfunUi = @(y, z) feval(fitfun, U(:, i), y, z);
		optionsU2i = options2;
		optionsU2i.initial = innerState{i};
		optionsU2i.initial.X = V2(:, :, i);
		optionsU2i.initial.f = [];
		optionsU2i.initial.A = []; % disabled external option archive
		optionsU2i.TolX = max(options2.TolX, innerTolX);
		
		for j = 1 : NP2
			optionsU2i.initial.innerState{j}.X = V3(:, :, j, i);
			optionsU2i.initial.innerState{j}.f = [];
		end
		
		[innerUbest1(:, i), innerUbest2(:, i), innerFbest, innerOutUi] = ...
			feval(innerSolver, innerFitfunUi, innerMaxfunevalsX, ...
			lb2, ub2, ...
			lb3, ub3, ...
			optionsU2i, options3);
		
		U_Converged_FEs(i) = innerOutUi.fes(end);
		counteval = counteval + innerOutUi.fes(end);
		fui = -innerFbest;
		
		% Replacement
		if fui < f(i)
			f(i) = fui;
			X(:, i) = U(:, i);
			innerXbest1(:, i) = innerUbest1(:, i);
			innerXbest2(:, i) = innerUbest2(:, i);
			innerState{i} = innerOutUi.final;
			successRate = successRate + 1 / NP1;
			S_F(Succ_Counter + 1) = F(i);
			S_CR(Succ_Counter + 1) = CR(i);
			Succ_Counter = Succ_Counter + 1;	
			FailedIteration = false;
		else
			innerState{i} = innerOutXi.final;
		end
	end
	
	% Display
	if isDisplayIter
		displayitermessages(...
			[X; innerXbest1], [U; innerUbest1], f, countiter, XX, YY, ZZ, ...
			'counteval', counteval, ...
			'successRate', successRate, ...
			'X_Converged_FEs', mean(X_Converged_FEs), ...
			'U_Converged_FEs', mean(U_Converged_FEs), ...
			'mu_F', mu_F, ...
			'mu_CR', mu_CR);
		
		display_inner_info(innerState);
	end
	
	% Update CR and F
	if Succ_Counter > 0
		mu_CR = (1-w) * mu_CR + w * mean(S_CR(1 : Succ_Counter));
		mu_F = (1-w) * mu_F + w * sum(S_F(1 : Succ_Counter).^2) / sum(S_F(1 : Succ_Counter));
	else
		mu_F = (1-w) * mu_F;
	end
	
	% Sort
	[f, fidx] = sort(f);
	X = X(:, fidx);
	innerState = innerState(fidx);
	
	% Record
	out = updateoutput(out, X, f, counteval, ...
		'innerFstd', computeInnerFstd(innerState), ...
		'innerMeanXstd', computeInnerMeanXstd(innerState), ...
		'successRate', successRate, ...
		'X_Converged_FEs', mean(X_Converged_FEs), ...
		'U_Converged_FEs', mean(U_Converged_FEs), ...
		'mu_F', mu_F, ...
		'mu_CR', mu_CR);
	
	% Iteration counter
	countiter = countiter + 1;
	
	% Stagnation iteration
	if FailedIteration
		countStagnation = countStagnation + 1;
	else
		countStagnation = 0;
	end
end

[fbest, fbestidx1] = min(f);
fbest = -fbest;
xbest1 = X(:, fbestidx1);
[~, fbestidx2] = min(innerState{fbestidx1}.f);
xbest2 = innerState{fbestidx1}.X(:, fbestidx2);
[~, fbestidx3] = min(innerState{fbestidx1}.innerState{fbestidx2}.f);
xbest3 = innerState{fbestidx1}.innerState{fbestidx2}.X(:, fbestidx3);

final.mu_F = mu_F;
final.mu_CR = mu_CR;
final.innerState = innerState;
out = finishoutput(out, X, f, counteval, 'final', final, ...
	'innerFstd', computeInnerFstd(innerState), ...
	'innerMeanXstd', computeInnerMeanXstd(innerState), ...
	'successRate', successRate, ...
	'X_Converged_FEs', mean(X_Converged_FEs), ...
	'U_Converged_FEs', mean(U_Converged_FEs), ...
	'mu_F', mu_F, ...
	'mu_CR', mu_CR);
end