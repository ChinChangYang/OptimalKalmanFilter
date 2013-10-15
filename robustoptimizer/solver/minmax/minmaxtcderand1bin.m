function [xbest1, xbest2, fbest, out] = minmaxtcderand1bin(fitfun, ...
	maxfunevals, lb1, ub1, lb2, ub2, options1, options2)
% MINMAXTCDERAND1BIN Min-Max Tracer Coevolutionary DE/rand/1
% MINMAXTCDERAND1BIN(fitfun, maxfunevals1, lb1, ub1, lb2, ub2) minimizes the
% function fitfun1 associated with a maximizer among box limitations [lb1,
% ub1] of minimizers and [lb2, ub2] of maximizers for the maximal function
% evaluations maxfunevals1.
% MINMAXTCDERAND1BIN(..., options1) minimizes the function with the given
% options Options1 for the 1st layer.
% MINMAXTCDERAND1BIN(..., options1, options2) minimize the function with the
% given options Options2 for the 2nd layer.
%
% Property:
% * Resumable solver (confirmed by test_resumable_minmax.m)
if nargin <= 6
	options1 = [];
end

if nargin <= 7
	options2 = [];
end

D1 = numel(lb1);
D2 = numel(lb2);

% Default options for Layer 1
defaultOptions1.dimensionFactor = 10;
defaultOptions1.F = 0.9;
defaultOptions1.CR = 0.99;
defaultOptions1.Display = 'off';
defaultOptions1.RecordPoint = 100;
defaultOptions1.ftarget = -Inf;
defaultOptions1.TolX = 0;
defaultOptions1.TolFun = 0;
defaultOptions1.TolStagnationIteration = 20;
defaultOptions1.TolX_DecayRate = 0.5;
defaultOptions1.InnerSolver = 'derand1bin';
defaultOptions1.initial.X = [];
defaultOptions1.initial.f = [];
defaultOptions1.initial.innerState = [];

defaultOptions1.TolCon = 1e-6;
defaultOptions1.nonlcon = [];
defaultOptions1.initial.cv = []; % Constraint violation measure
defaultOptions1.initial.nvc = []; % Number of violated constraints

options1 = setdefoptions(options1, defaultOptions1);

% Default options for Layer 2
defaultOptions2.dimensionFactor = 10;
defaultOptions2.F = 0.9;
defaultOptions2.Display = 'off';
defaultOptions2.RecordPoint = 0;
defaultOptions2.TolFun = 0;
defaultOptions2.TolX = 0;
options2 = setdefoptions(options2, defaultOptions2);

% Initialize algorithmic variables
dimensionFactor = max(1, options1.dimensionFactor);
F = options1.F;
CR = options1.CR;
isDisplayIter = strcmp(options1.Display, 'iter');
RecordPoint = max(0, floor(options1.RecordPoint));
TolFun = options1.TolFun;
TolX = options1.TolX;
TolStagnationIteration = options1.TolStagnationIteration;
TolX_DecayRate = options1.TolX_DecayRate;
innerSolver = options1.InnerSolver;
TolCon = options1.TolCon;
nonlcon = options1.nonlcon;

X = options1.initial.X;
f = options1.initial.f;
innerState = options1.initial.innerState;
existInnerState = ~isempty(innerState);
cv = options1.initial.cv;
nvc = options1.initial.nvc;

NP1 = ceil(dimensionFactor * D1);
NP2 = ceil(options2.dimensionFactor * D2);

% Initialize contour data
if isDisplayIter
	[XX, YY, ZZ] = minmaxcontourdata(D1, lb1, ub1, lb2, ub2, fitfun);
% 	load('minmaxcontourdata.mat');
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
innerXbest = zeros(D2, NP1);
innerUbest = innerXbest;
V = X;
U = X;
V2 = zeros(D2, NP2, NP1);
fu = zeros(1, NP1);
innerOutX = cell(1, NP1);
innerOutU = cell(1, NP1);
cv_u = zeros(1, NP1);
nvc_u = zeros(1, NP1);

out = initoutput(RecordPoint, D1, NP1, maxfunevals, ...
	'innerFstd', ...
	'innerMeanXstd', ...
	'successRate', ...
	'X_Converged_FEs', ...
	'U_Converged_FEs');

% Constraint violation measure
if isempty(cv) || isempty(nvc)
	cv = zeros(1, NP1);
	nvc = zeros(1, NP1);
	
	if ~isempty(nonlcon)		
		for i = 1 : NP1
			[c, ceq] = feval(nonlcon, X(:, i));
			cv(i) = sum(c) + sum(ceq);
			nvc(i) = sum(c > 0) + sum(ceq > 0);
		end
	end
end

% Evaluation
if isempty(f)
	f = zeros(1, NP1);
	innerMaxfunevalsX = NP2;
	
	for i = 1 : NP1
		if nvc(i) > 0
			f(i) = inf;
		else
			innerFitfun = @(y) -feval(fitfun, X(:, i), y);
			optionsX2i = options2;
			
			if existInnerState
				optionsX2i.initial = innerState{i};
			end
			
			[innerXbest(:, i), innerFbest, innerOut] = ...
				feval(innerSolver, innerFitfun, ...
				lb2, ub2, innerMaxfunevalsX, optionsX2i);
			
			counteval = counteval + innerOut.fes(end);
			f(i) = -innerFbest;
			innerState{i} = innerOut.final;
		end
	end
end

% Sort
pf = zeros(1, NP1);
nf = f;
nf(isinf(nf)) = [];
nfmax = max(nf);
nfmin = min(nf);
ncv = cv;
ncvmax = max(ncv);
ncvmin = min(ncv);

for i = 1 : NP1
	if nvc(i) == 0
		pf(i) = (f(i) - nfmin) / (nfmax - nfmin + eps);
	else
		pf(i) = nvc(i) + (ncv(i) - ncvmin) / (ncvmax - ncvmin + eps);
	end
end

[~, pfidx] = sort(pf);
f = f(pfidx);
X = X(:, pfidx);
innerState = innerState(pfidx);
cv = cv(pfidx);
nvc = nvc(pfidx);

% Display
if isDisplayIter
	displayitermessages([X; innerXbest], [X; innerXbest], f, countiter, ...
		XX, YY, ZZ, 'counteval', counteval, ...
		'successRate', successRate);
	
	display_inner_info(innerState);
end

% Record minimal function values
out = updateoutput(out, X, f, counteval, ...
	'innerFstd', computeInnerFstd(innerState), ...
	'innerMeanXstd', computeInnerMeanXstd(innerState), ...
	'successRate', successRate, ...
	'X_Converged_FEs', mean(X_Converged_FEs), ...
	'U_Converged_FEs', mean(U_Converged_FEs));

countiter = countiter + 1;

while true
	% Termination conditions
	outofmaxfunevals = counteval >= maxfunevals;
	solutionconvergence = isConverged(X, TolX);
	stagnation = countStagnation >= TolStagnationIteration;
	fitnessconvergence = isConverged(f, TolFun) && ...
		isConverged(cv, TolCon);
	
	% Convergence conditions
	if outofmaxfunevals || fitnessconvergence || solutionconvergence ...
			|| stagnation
		break;
	end
	
	% Mutation
	for i = 1 : NP1
		% Try generating V within bounds
		for retry_within_bounds = 1 : NP1
			
			% Generate r1
			r1 = floor(1 + NP1 * rand);
			
			% Generate r2
			for retry = 1 : 3
				r2 = floor(1 + NP1 * rand);
				if ~(all(X(:, r1) == X(:, r2)))
					break;
				end
			end
			
			% Generate r3
			for retry = 1 : 3
				r3 = floor(1 + NP1 * rand);
				if ~(all(X(:, r1) == X(:, r3)) || all(X(:, r2) == X(:, r3)))
					break;
				end
			end
			
			% Generate Fi
			Fi = F + 0.01 * randn;
			
			% Generate Vi
			V(:, i) = X(:, r1) + Fi .* (X(:, r2) - X(:, r3));
			
			% Check boundary
			if all(V(:, i) >= lb1) && all(V(:, i) <= ub1)
				
				if ~isempty(innerState{r1}) && ...
						~isempty(innerState{r2}) && ...
						~isempty(innerState{r3})
					
					% Prediction
					V2(:, :, i) = innerState{r1}.X + ...
						Fi .* (innerState{r2}.X - innerState{r3}.X);
					
					% Perturbation
					V2(:, :, i) = V2(:, :, i) .* (1 + 100 * eps * randn(D2, NP2));
				else
					for j = 1 : NP2
						V2(:, j, i) = lb2 + (ub2 - lb2) .* rand(D2, 1);
					end
				end
				
				break;
			end
		end
	end
	
	% Crossover
	for i = 1 : NP1
		jrand = floor(1 + D1 * rand);
		
		for j = 1 : D1
			if rand < CR || j == jrand
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
	end
	
	% Constraint violation measure
	if ~isempty(nonlcon)
		for i = 1 : NP1
			[c, ceq] = feval(nonlcon, U(:, i));
			cv_u(i) = sum(c) + sum(ceq);
			nvc_u(i) = sum(c > 0) + sum(ceq > 0);
		end
	end
	
	% TolX Decision
	innerTolX = TolX_DecayRate * computeInnerMeanXstd(innerState);
		
	% Selection
	innerMaxfunevalsX = 20 * NP2;
	for i = 1 : NP1
		% Compute fxi, f(i)
		if nvc(i) == 0
			innerFitfunXi = @(X2) -feval(fitfun, X(:, i), X2);
			optionsX2i = options2;
			optionsX2i.initial = innerState{i};
			optionsX2i.TolX = max(options2.TolX, innerTolX);
			
			[innerXbest(:, i), innerFbest, innerOutX{i}] = ...
				feval(innerSolver, innerFitfunXi, ...
				lb2, ub2, ...
				innerMaxfunevalsX, optionsX2i);
			
			X_Converged_FEs(i) = innerOutX{i}.fes(end);
			counteval = counteval + innerOutX{i}.fes(end);
			f(i) = -innerFbest;
		else
			innerOutX{i}.final = [];
		end
		
		% Compute fui
		if nvc_u(i) == 0
			fitfunU2i = @(U2) -feval(fitfun, U(:, i), U2);
			optionsU2i = options2;
			optionsU2i.initial = innerState{i};
			optionsU2i.initial.X = V2(:, :, i);
			optionsU2i.initial.f = [];
			optionsU2i.initial.cv = [];
			optionsU2i.initial.nvc = [];
			optionsU2i.TolX = max(options2.TolX, innerTolX);
			
			[innerUbest(:, i), innerFbest, innerOutU{i}] = ...
				feval(innerSolver, fitfunU2i, ...
				lb2, ub2, ...
				innerMaxfunevalsX, optionsU2i);
			
			U_Converged_FEs(i) = innerOutU{i}.fes(end);
			counteval = counteval + innerOutU{i}.fes(end);
			fu(i) = -innerFbest;
		else
			fu(i) = inf;
			innerOutU{i}.final = [];
		end
	end
	
	% Replacement
	successRate = 0;
	FailedIteration = true;
	for i = 1 : NP1		
		if nvc(i) == 0 && nvc_u(i) == 0
			if fu(i) < f(i)
				u_selected = true;
			else
				u_selected = false;
			end
		elseif nvc(i) > nvc_u(i)
			u_selected = true;
		elseif nvc(i) < nvc_u(i)
			u_selected = false;
		else % nvc(i) == nvc_u(i) && nvc(i) ~= 0 && nvc_u(i) ~= 0
			if cv(i) > cv_u(i)
				u_selected = true;
			else
				u_selected = false;
			end
		end
		
		if u_selected
			cv(i) = cv_u(i);
			nvc(i) = nvc_u(i);
			f(i) = fu(i);
			X(:, i) = U(:, i);
			innerXbest(:, i) = innerUbest(:, i);
			innerState{i} = innerOutU{i}.final;
			successRate = successRate + 1 / NP1;
			FailedIteration = false;
		else
			innerState{i} = innerOutX{i}.final;
		end
	end
	
	% Display
	if isDisplayIter
		displayitermessages([X; innerXbest], [U; innerUbest], f, countiter, ...
			XX, YY, ZZ, 'counteval', counteval, ...
			'successRate', successRate, ...
			'X_Converged_FEs', mean(X_Converged_FEs), ...
			'U_Converged_FEs', mean(U_Converged_FEs));
		
		display_inner_info(innerState);
	end
	
	% Sort
	nf = f;
	nf(isinf(nf)) = [];
	nfmax = max(nf);
	nfmin = min(nf);
	ncv = cv;
	ncvmax = max(ncv);
	ncvmin = min(ncv);
	
	for i = 1 : NP1
		if nvc(i) == 0
			pf(i) = (f(i) - nfmin) / (nfmax - nfmin + eps);
		else
			pf(i) = nvc(i) + (ncv(i) - ncvmin) / (ncvmax - ncvmin + eps);
		end
	end
	
	[~, pfidx] = sort(pf);
	f = f(pfidx);
	X = X(:, pfidx);
	innerState = innerState(pfidx);
	cv = cv(pfidx);
	nvc = nvc(pfidx);
	
	% Record
	out = updateoutput(out, X, f, counteval, ...
		'innerFstd', computeInnerFstd(innerState), ...
		'innerMeanXstd', computeInnerMeanXstd(innerState), ...
		'successRate', successRate, ...
		'X_Converged_FEs', mean(X_Converged_FEs), ...
		'U_Converged_FEs', mean(U_Converged_FEs));
	
	% Iteration counter
	countiter = countiter + 1;
	
	% Stagnation iteration
	if FailedIteration
		countStagnation = countStagnation + 1;
	else
		countStagnation = 0;
	end
end

[fbest, fbestidx] = min(f);
xbest1 = X(:, fbestidx);
[~, fbestidx2] = min(innerState{fbestidx}.f);
xbest2 = innerState{fbestidx}.X(:, fbestidx2);

final.innerState = innerState;
final.cv = cv;
final.nvc = nvc;

out = finishoutput(out, X, f, counteval, 'final', final, ...
	'innerFstd', computeInnerFstd(innerState), ...
	'innerMeanXstd', computeInnerMeanXstd(innerState), ...
	'successRate', successRate, ...
	'X_Converged_FEs', mean(X_Converged_FEs), ...
	'U_Converged_FEs', mean(U_Converged_FEs));
end