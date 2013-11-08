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
defaultOptions1.CR = 0.5;
defaultOptions1.delta_F = 0.1;
defaultOptions1.delta_CR = 0.1;
defaultOptions1.p = 0.05;
defaultOptions1.w = 0.1;
defaultOptions1.Display = 'off';
defaultOptions1.RecordPoint = 100;
defaultOptions1.TolX = 0;
defaultOptions1.TolFun = 0;
defaultOptions1.TolStagnationIteration = 20;
defaultOptions1.InnerSolver = 'minmaxtcjadebin';
defaultOptions1.initial.X = [];
defaultOptions1.initial.f = [];
defaultOptions1.initial.mu_F = [];
defaultOptions1.initial.mu_CR = [];
defaultOptions1.initial.innerState = [];

defaultOptions1.TolCon = 1e-6;
defaultOptions1.nonlcon = [];
defaultOptions1.innerMaxIter = 100;
defaultOptions1.reinitFactor = 0.1;
defaultOptions1.migrateFactor = 0.8;
defaultOptions1.archiveSizeFactor = 4;

options1 = setdefoptions(options1, defaultOptions1);

% Default options for Layer 2
defaultOptions2.dimensionFactor = 10;
defaultOptions2.Display = 'off';
defaultOptions2.RecordPoint = 0;
defaultOptions2.TolFun = 0;
defaultOptions2.TolX = 0;
defaultOptions2.InnerSolver = 'jadebin';
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
innerSolver = options1.InnerSolver;
TolCon = options1.TolCon;
nonlcon = options1.nonlcon;
innerMaxIter = options1.innerMaxIter;
reinitFactor = options1.reinitFactor;
migrateFactor = options1.migrateFactor;
archiveSizeFactor = options1.archiveSizeFactor;

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
innerXbest1 = zeros(D2, NP1);
innerXbest2 = zeros(D3, NP1);
innerUbest1 = innerXbest1;
innerUbest2 = innerXbest2;
V = X;
U = X;
V2 = zeros(D2, NP2, NP1);
V3 = zeros(D3, NP3, NP2, NP1);
pbest_size = p * NP1;
fu = zeros(1, NP1);
innerOutX = cell(1, NP1);
innerOutU = cell(1, NP1);
cm_u = zeros(1, NP1);
nc_u = zeros(1, NP1);
reinitNP2 = ceil(reinitFactor * NP2);
reinitNP3 = ceil(reinitFactor * NP3);
archiveSize = ceil(archiveSizeFactor * NP1);

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
	innerMaxfunevalsX = innerMaxIter * NP2 * NP3;
	
	parfor i = 1 : NP1
		innerFitfun = @(y, z) feval(fitfun, X(:, i), y, z);	
		optionsX2i = options2;
		
		if ~isempty(nonlcon)
			innerNonlcon = @(y, z) feval(nonlcon, X(:, i), y, z);
			optionsX2i.nonlcon = innerNonlcon;
		end
		
		if existInnerState
			optionsX2i.initial = innerState{i};
		end
		
		[innerXbest1(:, i), innerXbest2(:, i), innerFbest, innerOut] = ...
			feval(innerSolver, innerFitfun, innerMaxfunevalsX, lb2, ub2, ...
			lb3, ub3, optionsX2i, options3);
		
		counteval = counteval + innerOut.fes(end);
		f(i) = -innerFbest;
		innerState{i} = innerOut.final;
	end
end

% Initialize archive
archive1 = repmat(innerXbest1, 1, archiveSizeFactor);
archive2 = repmat(innerXbest2, 1, archiveSizeFactor);

% Constraint violation measure
cm = zeros(1, NP1);
nc = zeros(1, NP1);

for i = 1 : NP1
	clb = lb1 - X(:, i);
	cub = X(:, i) - ub1;
	cm(i) = sum(clb(clb > 0)) + sum(cub(cub > 0));
	nc(i) = sum(clb > 0) + sum(cub > 0);
end

for i = 1 : NP1
	clb = lb2 - innerXbest1(:, i);
	cub = innerXbest1(:, i) - ub2;
	cm(i) = cm(i) + sum(clb(clb > 0)) + sum(cub(cub > 0));
	nc(i) = nc(i) + sum(clb > 0) + sum(cub > 0);
end

for i = 1 : NP1
	clb = lb3 - innerXbest2(:, i);
	cub = innerXbest2(:, i) - ub3;
	cm(i) = cm(i) + sum(clb(clb > 0)) + sum(cub(cub > 0));
	nc(i) = nc(i) + sum(clb > 0) + sum(cub > 0);
end

if ~isempty(nonlcon)
	for i = 1 : NP1
		[cx, ceqx] = ...
			feval(nonlcon, X(:, i), innerXbest1(:, i), innerXbest2(:, i));
		
		cm(i) = cm(i) + sum(cx(cx > 0)) + sum(ceqx(ceqx > 0));
		nc(i) = nc(i) + sum(cx > 0) + sum(ceqx > 0);
	end
end

% Sort
pf = zeros(1, NP1);
nf = f;
nf(isinf(nf)) = [];
nfmax = max(nf);
nfmin = min(nf);
cmmax = max(cm);
cmmin = min(cm);

for i = 1 : NP1
	if nc(i) == 0
		pf(i) = (f(i) - nfmin) / (nfmax - nfmin + eps);
	else
		pf(i) = nc(i) + (cm(i) - cmmin) / (cmmax - cmmin + eps);
	end
end

[pf, pfidx] = sort(pf);
f = f(pfidx);
X = X(:, pfidx);
innerXbest1 = innerXbest1(:, pfidx);
innerXbest2 = innerXbest2(:, pfidx);
innerState = innerState(pfidx);
cm = cm(pfidx);
nc = nc(pfidx);

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
	if all(isinf(f))
		displayitermessages([X; innerXbest1], [X; innerXbest1], ...
			cm, countiter, ...
			XX, YY, ZZ, 'counteval', counteval, ...
			'successRate', successRate, ...
			'mu_F', mu_F, ...
			'mu_CR', mu_CR);
	else
		displayitermessages([X; innerXbest1], [X; innerXbest1], ...
			f(~isinf(f)), countiter, ...
			XX, YY, ZZ, 'counteval', counteval, ...
			'successRate', successRate, ...
			'mu_F', mu_F, ...
			'mu_CR', mu_CR);
	end
	
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
	fitnessconvergence = isConverged(f, TolFun) && isConverged(cm, TolCon);
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
	
	% Prediction
	anyEmptyInnerState = false;
	for i = 1 : NP1
		if isempty(innerState{i})
			anyEmptyInnerState = true;
			break;
		end
		
		for j = 1 : NP2
			if isempty(innerState{i}.innerState{j})
				anyEmptyInnerState = true;
				break;
			end
		end
		
		if anyEmptyInnerState
			break;
		end
	end
	
	if ~anyEmptyInnerState		
		for i = 1 : NP1			
			% Copy from itselfs individuals
			for j = 1 : NP2
				V2(:, j, i) = innerState{i}.X(:, j);
				
				for k = 1 : NP3
					V3(:, k, j, i) = innerState{i}.innerState{j}.X(:, k);
				end
			end
			
			% Copy from archive
			migrationNP2 = ceil(migrateFactor * NP2);
			beginIndex1 = NP2 - migrationNP2 - reinitNP2;
			for j = 1 : migrationNP2
				r = floor(archiveSize * rand + 1);
				V2(:, beginIndex1 + j, i) = archive1(:, r);
			end
			
			migrationNP3 = ceil(migrateFactor * NP3);
			beginIndex2 = NP3 - migrationNP3 - reinitNP3;
			for j = 1 : NP2
				for k = 1 : migrationNP3
					r = floor(archiveSize * rand + 1);
					V3(:, beginIndex2 + k, j, i) = archive2(:, r);
				end
			end
			
			% Reinitialization
			for j = (NP2 - reinitNP2 + 1) : NP2
				V2(:, j, i) = ...
					lb2 + (ub2 - lb2) .* rand(D2, 1);
			end
			
			for j = 1 : NP2
				for k = (NP3 - reinitNP3 + 1) : NP3
					V3(:, k, j, i) = ...
						lb3 + (ub3 - lb3) .* rand(D3, 1);
				end
			end
		end
	else		
		for i = 1 : NP1
			for j = 1 : NP2
				V2(:, j, i) = lb2 + (ub2 - lb2) .* rand(D2, 1);
				
				for k = 1 : NP3
					V3(:, k, j, i) = lb3 + (ub3 - lb3) .* rand(D3, 1);
				end
			end
		end
	end
	
	% Selection
	innerMaxfunevalsX = innerMaxIter * NP2 * NP3;
	parfor i = 1 : NP1
		% Compute fxi, f(i)
		innerFitfunXi = @(y, z) feval(fitfun, X(:, i), y, z);
		optionsX2i = options2;
		optionsX2i.initial = innerState{i};		
		optionsX2i.initial.X = V2(:, :, i);
		optionsX2i.initial.f = [];
		optionsX2i.initial.mu_F = [];
		optionsX2i.initial.mu_CR = [];
		optionsX2i.initial.cm = [];
		optionsX2i.initial.nc = [];
		
		for j = 1 : NP2
			optionsX2i.initial.innerState{j} = [];
			optionsX2i.initial.innerState{j}.X = V3(:, :, j, i);
			optionsX2i.initial.innerState{j}.f = [];
			optionsX2i.initial.innerState{j}.A = [];
			optionsX2i.initial.innerState{j}.mu_F = [];
			optionsX2i.initial.innerState{j}.mu_CR = [];
			optionsX2i.initial.innerState{j}.cm = [];
			optionsX2i.initial.innerState{j}.nc = [];
		end
		
		if ~isempty(nonlcon)
			optionsX2i.nonlcon = @(y, z) feval(nonlcon, X(:, i), y, z);
		end
		
		[innerXbest1(:, i), innerXbest2(:, i), innerFbest, innerOutX{i}] = ...
			feval(innerSolver, innerFitfunXi, innerMaxfunevalsX, ...
			lb2, ub2, ...
			lb3, ub3, ...
			optionsX2i, options3);
		
		X_Converged_FEs(i) = innerOutX{i}.fes(end);
		counteval = counteval + innerOutX{i}.fes(end);
		f(i) = -innerFbest;
		
		% Compute fui
		innerFitfunUi = @(y, z) feval(fitfun, U(:, i), y, z);
		optionsU2i = options2;
		optionsU2i.initial = innerState{i};
		optionsU2i.initial.X = V2(:, :, i);
		optionsU2i.initial.f = [];
		optionsU2i.initial.A = [];
		optionsU2i.initial.mu_F = [];
		optionsU2i.initial.mu_CR = [];
		optionsU2i.initial.cm = [];
		optionsU2i.initial.nc = [];
		
		for j = 1 : NP2
			optionsU2i.initial.innerState{j} = [];
			optionsU2i.initial.innerState{j}.X = V3(:, :, j, i);
			optionsU2i.initial.innerState{j}.f = [];
			optionsU2i.initial.innerState{j}.A = [];
		end
		
		[innerUbest1(:, i), innerUbest2(:, i), innerFbest, innerOutU{i}] = ...
			feval(innerSolver, innerFitfunUi, innerMaxfunevalsX, ...
			lb2, ub2, ...
			lb3, ub3, ...
			optionsU2i, options3);
		
		U_Converged_FEs(i) = innerOutU{i}.fes(end);
		counteval = counteval + innerOutU{i}.fes(end);
		fu(i) = -innerFbest;
	end
	
	% Constraint violation measure
	for i = 1 : NP1
		clb = lb1 - X(:, i);
		cub = X(:, i) - ub1;
		cm(i) = sum(clb(clb > 0)) + sum(cub(cub > 0));
		nc(i) = sum(clb > 0) + sum(cub > 0);
		
		clb = lb1 - U(:, i);
		cub = U(:, i) - ub1;
		cm_u(i) = sum(clb(clb > 0)) + sum(cub(cub > 0));
		nc_u(i) = sum(clb > 0) + sum(cub > 0);
	end
	
	for i = 1 : NP1
		clb = lb2 - innerXbest1(:, i);
		cub = innerXbest1(:, i) - ub2;
		cm(i) = cm(i) + sum(clb(clb > 0)) + sum(cub(cub > 0));
		nc(i) = nc(i) + sum(clb > 0) + sum(cub > 0);
				
		clb = lb2 - innerUbest1(:, i);
		cub = innerUbest1(:, i) - ub2;
		cm_u(i) = cm_u(i) + sum(clb(clb > 0)) + sum(cub(cub > 0));
		nc_u(i) = nc_u(i) + sum(clb > 0) + sum(cub > 0);
	end
	
	for i = 1 : NP1
		clb = lb3 - innerXbest2(:, i);
		cub = innerXbest2(:, i) - ub3;
		cm(i) = cm(i) + sum(clb(clb > 0)) + sum(cub(cub > 0));
		nc(i) = nc(i) + sum(clb > 0) + sum(cub > 0);
				
		clb = lb3 - innerUbest2(:, i);
		cub = innerUbest2(:, i) - ub3;
		cm_u(i) = cm_u(i) + sum(clb(clb > 0)) + sum(cub(cub > 0));
		nc_u(i) = nc_u(i) + sum(clb > 0) + sum(cub > 0);
	end
	
	if ~isempty(nonlcon)
		for i = 1 : NP1
			[cx, ceqx] = feval(...
				nonlcon, ...
				X(:, i), ...
				innerXbest1(:, i), ...
				innerXbest2(:, i));
			
			cm(i) = cm(i) + sum(cx(cx > 0)) + sum(ceqx(ceqx > 0));
			nc(i) = nc(i) + sum(cx > 0) + sum(ceqx > 0);
			
			[cu, cequ] = feval(...
				nonlcon, ...
				U(:, i), ...
				innerUbest1(:, i), ...
				innerUbest2(:, i));
			
			cm_u(i) = cm_u(i) + sum(cu(cu > 0)) + sum(cequ(cequ > 0));
			nc_u(i) = nc_u(i) + sum(cu > 0) + sum(cequ > 0);
		end
	end
		
	% Replacement
	successRate = 0;
	FailedIteration = true;
	for i = 1 : NP1
		if nc(i) == 0 && nc_u(i) == 0
			if fu(i) < f(i)
				u_selected = true;
			else
				u_selected = false;
			end
		elseif nc(i) > nc_u(i)
			u_selected = true;
		elseif nc(i) < nc_u(i)
			u_selected = false;
		else % nvc(i) == nvc_u(i) && nvc(i) ~= 0 && nvc_u(i) ~= 0
			if cm(i) > cm_u(i)
				u_selected = true;
			else
				u_selected = false;
			end
		end
		
		if u_selected
			cm(i) = cm_u(i);
			nc(i) = nc_u(i);
			f(i) = fu(i);
			X(:, i) = U(:, i);
			innerXbest1(:, i) = innerUbest1(:, i);
			innerXbest2(:, i) = innerUbest2(:, i);
			innerState{i} = innerOutU{i}.final;
			successRate = successRate + 1 / NP1;
			S_F(Succ_Counter + 1) = F(i);
			S_CR(Succ_Counter + 1) = CR(i);
			Succ_Counter = Succ_Counter + 1;	
			FailedIteration = false;
		else
			innerState{i} = innerOutX{i}.final;
		end
	end
	
	% Display
	if isDisplayIter
		if all(isinf(f))
			displayitermessages([X; innerXbest1], [X; innerXbest1], ...
				cm, countiter, ...
				XX, YY, ZZ, 'counteval', counteval, ...
				'successRate', successRate, ...
				'mu_F', mu_F, ...
				'mu_CR', mu_CR);
		else
			displayitermessages([X; innerXbest1], [X; innerXbest1], ...
				f(~isinf(f)), countiter, ...
				XX, YY, ZZ, 'counteval', counteval, ...
				'successRate', successRate, ...
				'mu_F', mu_F, ...
				'mu_CR', mu_CR);
		end
		
		display_inner_info(innerState);
	end
	
	% Update archive
	randindex = randperm(archiveSize);
	archive1(:, randindex(1 : NP1)) = innerXbest1;
	archive2(:, randindex(1 : NP1)) = innerXbest2;
	
	% Update CR and F
	if Succ_Counter > 0
		mu_CR = (1-w) * mu_CR + w * mean(S_CR(1 : Succ_Counter));
		mu_F = (1-w) * mu_F + w * sum(S_F(1 : Succ_Counter).^2) / sum(S_F(1 : Succ_Counter));
	else
		mu_F = (1-w) * mu_F;
	end
	
	% Sort
	nf = f;
	nf(isinf(nf)) = [];
	nfmax = max(nf);
	nfmin = min(nf);
	cmmax = max(cm);
	cmmin = min(cm);
	
	for i = 1 : NP1
		if nc(i) == 0
			pf(i) = (f(i) - nfmin) / (nfmax - nfmin + eps);
		else
			pf(i) = nc(i) + (cm(i) - cmmin) / (cmmax - cmmin + eps);
		end
	end
	
	[pf, pfidx] = sort(pf);
	f = f(pfidx);
	X = X(:, pfidx);
	innerXbest1 = innerXbest1(:, pfidx);
	innerXbest2 = innerXbest2(:, pfidx);
	innerState = innerState(pfidx);
	cm = cm(pfidx);
	nc = nc(pfidx);
	
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

xbest1 = X(:, 1);
xbest2 = innerState{1}.X(:, 1);
xbest3 = innerState{1}.innerState{1}.X(:, 1);
fbest = -f(1);

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