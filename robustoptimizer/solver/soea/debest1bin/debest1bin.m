function [xmin, fmin, out] = debest1bin(fitfun, lb, ub, maxfunevals, options)
% DEBEST1BIN Classical DE/best/1/bin
% DEBEST1BIN(fitfun, lb, ub, maxfunevals) minimize the function fitfun in
% box constraints [lb, ub] with the maximal function evaluations
% maxfunevals.
% DEBEST1BIN(..., options) minimize the function by solver options.
if nargin <= 4
	options = [];
end

defaultOptions.dimensionFactor = 5;
defaultOptions.CR = 0.5;
defaultOptions.F = 0.7;
defaultOptions.Display = 'off';
defaultOptions.RecordPoint = 100;
defaultOptions.ftarget = -Inf;
defaultOptions.TolFun = eps;
defaultOptions.TolX = 100 * eps;
defaultOptions.TolStagnationIteration = 20;
defaultOptions.initial.X = [];
defaultOptions.initial.f = [];

defaultOptions.TolCon = 1e-6;
defaultOptions.nonlcon = [];
defaultOptions.initial.cv = []; % Constraint violation measure
defaultOptions.initial.nvc = []; % Number of violated constraints

options = setdefoptions(options, defaultOptions);
dimensionFactor = options.dimensionFactor;
CR = options.CR;
F = options.F;
isDisplayIter = strcmp(options.Display, 'iter');
RecordPoint = max(1, floor(options.RecordPoint));
ftarget = options.ftarget;
TolFun = options.TolFun;
TolX = options.TolX;
TolStagnationIteration = options.TolStagnationIteration;
TolCon = options.TolCon;
nonlcon = options.nonlcon;

if ~isempty(options.initial)
	X = options.initial.X;
	f = options.initial.f;
	cv = options.initial.cv;
	nvc = options.initial.nvc;
else
	X = [];
	f = [];
	cv = [];
	nvc = [];
end

D = numel(lb);
if isempty(X)
	NP = ceil(dimensionFactor * D);
else
	[~, NP] = size(X);
end

% Initialize variables
counteval = 0;
countiter = 1;
countStagnation = 0;
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

% Constraint violation measure
if isempty(cv) || isempty(nvc)
	cv = zeros(1, NP);
	nvc = zeros(1, NP);
	
	if ~isempty(nonlcon)		
		for i = 1 : NP
			[c, ceq] = feval(nonlcon, X(:, i));
			cv(i) = sum(c) + sum(ceq);
			nvc(i) = sum(c > 0) + sum(ceq > 0);
		end
	end
end

% Evaluation
if isempty(f)
	f = zeros(1, NP);
	for i = 1 : NP
		if nvc(i) > 0
			f(i) = inf;
		else
			f(i) = feval(fitfun, X(:, i));
			counteval = counteval + 1;
		end
	end
end

% Initialize variables
V = X;
U = X;
cv_u = cv;
nvc_u = nvc;

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
	outofmaxfunevals = counteval > maxfunevals - NP;
	reachftarget = min(f) <= ftarget;
	fitnessconvergence = isConverged(f, TolFun) && isConverged(cv, TolCon);
	solutionconvergence = isConverged(X, TolX);
	stagnation = countStagnation >= TolStagnationIteration;
	
	% Convergence conditions	
	if outofmaxfunevals || reachftarget || fitnessconvergence || ...
			solutionconvergence || stagnation
		break;
	end
	
	% Mutation
	[~, gbest] = min(f);
	for i = 1 : NP
		r1 = floor(1 + NP * rand);
		r2 = floor(1 + NP * rand);
		
		while r1 == r2
			r2 = floor(1 + NP * rand);
		end
		
		V(:, i) = X(:, gbest) + (F + 0.01 * randn) * (X(:, r1) - X(:, r2));
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
	
	% Constraint violation measure
	if ~isempty(nonlcon)
		for i = 1 : NP
			[c, ceq] = feval(nonlcon, U(:, i));
			cv_u(i) = sum(c) + sum(ceq);
			nvc_u(i) = sum(c > 0) + sum(ceq > 0);
		end
	end
	
	% Selection
	FailedIteration = true;
	for i = 1 : NP
		fui = inf;
		
		if nvc(i) == 0 && nvc_u(i) == 0
			fui = feval(fitfun, U(:, i));
			counteval = counteval + 1;
			
			if fui < f(i)
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
			f(i) = fui;
			X(:, i) = U(:, i);
			FailedIteration = false;
		end
	end
	
	% Record
	out = updateoutput(out, X, f, counteval);
	
	% Iteration counter
	countiter = countiter + 1;
	
	% Stagnation iteration
	if FailedIteration
		countStagnation = countStagnation + 1;
	else
		countStagnation = 0;
	end	
end

[fmin, fminidx] = min(f);
xmin = X(:, fminidx);

if fmin < out.bestever.fmin
	out.bestever.fmin = fmin;
	out.bestever.xmin = xmin;
end

final.cv = cv;
final.nvc = nvc;

out = finishoutput(out, X, f, counteval, 'final', final);
end