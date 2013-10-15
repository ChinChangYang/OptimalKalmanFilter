function [xmin, fmin, out] = jadebinwoa(fitfun, lb, ub, maxfunevals, options)
% JADEBINWOA JADE algorithm w/o external option archive
% JADEBINWOA(fitfun, lb, ub, maxfunevals) minimize the function fitfun in
% box constraints [lb, ub] with the maximal function evaluations
% maxfunevals.
% JADEBINWOA(..., options) minimize the function by solver options.
if nargin <= 4
	options = [];
end

defaultOptions.dimensionFactor = 5;
defaultOptions.F = 1.0;
defaultOptions.CR = 0.5;
defaultOptions.delta_CR = 0.1;
defaultOptions.delta_F = 0.1;
defaultOptions.p = 0.05;
defaultOptions.w = 0.1;
defaultOptions.Display = 'off';
defaultOptions.RecordPoint = 100;
defaultOptions.ftarget = -Inf;
defaultOptions.TolFun = eps;
defaultOptions.TolX = 100 * eps;
defaultOptions.TolStagnationIteration = 20;
defaultOptions.initial.X = [];
defaultOptions.initial.f = [];
defaultOptions.initial.mu_F = [];
defaultOptions.initial.mu_CR = [];

defaultOptions.TolCon = 1e-6;
defaultOptions.nonlcon = [];
defaultOptions.initial.cv = []; % Constraint violation measure
defaultOptions.initial.nvc = []; % Number of violated constraints

options = setdefoptions(options, defaultOptions);
dimensionFactor = max(1, options.dimensionFactor);
delta_CR = options.delta_CR;
delta_F = options.delta_F;
p = options.p;
w = options.w;
isDisplayIter = strcmp(options.Display, 'iter');
RecordPoint = max(0, floor(options.RecordPoint));
ftarget = options.ftarget;
TolFun = options.TolFun;
TolX = options.TolX;
TolStagnationIteration = options.TolStagnationIteration;
TolCon = options.TolCon;
nonlcon = options.nonlcon;

if ~isempty(options.initial)
	X = options.initial.X;
	f = options.initial.f;
	mu_CR = options.initial.mu_CR;
	mu_F = options.initial.mu_F;
	cv = options.initial.cv;
	nvc = options.initial.nvc;
else
	X = [];
	f = [];
	mu_CR = [];
	mu_F = [];
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

% Sort
[f, fidx] = sort(f);
X = X(:, fidx);

% mu_F
if isempty(mu_F)
	mu_F = options.F;
end

% mu_CR
if isempty(mu_CR)
	mu_CR = options.CR;
end

% Initialize variables
V = X;
U = X;
pbest_size = p * NP;
cv_u = cv;
nvc_u = nvc;

% Display
if isDisplayIter
	displayitermessages(...
		X, U, f, countiter, XX, YY, ZZ, 'mu_F', mu_F, 'mu_CR', mu_CR);
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
	
	% Scaling factor and crossover rate
	S_F = zeros(1, NP);
	S_CR = zeros(1, NP);
	CR = mu_CR + delta_CR * randn(1, NP);
	CR(CR > 1) = 1;
	CR(CR < 0) = 0;
	F = cauchyrnd(mu_F, delta_F, NP, 1);
	F(F > 1) = 1;
	
	for retry = 1 : 3
		if all(F > 0)
			break;
		end
		
		F(F <= 0) = cauchyrnd(mu_F, delta_F, sum(F <= 0), 1);
		F(F > 1) = 1;
	end
	
	F(F <= 0) = eps;
	
	A_Counter = 0;
	
	% Mutation
	for i = 1 : NP		
		for checkRetry = 1 : 3			
			% Generate pbest_idx
			for retry = 1 : 3
				pbest_idx = max(1, ceil(rand * pbest_size));
				if ~all(X(:, pbest_idx) == X(:, i))
					break;
				end
			end
			
			% Generate r1
			for retry = 1 : NP
				r1 = floor(1 + NP * rand);
				if i ~= r1
					break;
				end
			end
			
			% Generate r2
			for retry = 1 : NP * NP
				r2 = floor(1 + NP * rand);
				if ~(all(X(:, i) == X(:, r2)) || all(X(:, r1) == X(:, r2)))
					break;
				end
			end
							
			V(:, i) = X(:, i) + F(i) * (X(:, pbest_idx) - X(:, i) + X(:, r1) - X(:, r2));
			
			% Check boundary
			if all(V(:, i) > lb) && all(V(:, i) < ub)
				break;
			end
		end
	end
	
	for i = 1 : NP
		% Binominal Crossover
		jrand = floor(1 + D * rand);
		for j = 1 : D
			if rand < CR(i) || j == jrand
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
		displayitermessages(...
			X, U, f, countiter, XX, YY, ZZ, ...
			'mu_F', mu_F, 'mu_CR', mu_CR);
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
			S_CR(A_Counter + 1) = CR(i);
			S_F(A_Counter + 1) = F(i);
			A_Counter = A_Counter + 1;
			FailedIteration = false;
		end
	end
	
	% Update CR and F
	if A_Counter > 0
		mu_CR = (1-w) * mu_CR + w * mean(S_CR(1 : A_Counter));
		mu_F = (1-w) * mu_F + w * sum(S_F(1 : A_Counter).^2) / sum(S_F(1 : A_Counter));
	else
		mu_F = (1-w) * mu_F;
	end
	
	% Sort
	[f, fidx] = sort(f);
	X = X(:, fidx);
	
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

final.mu_F = mu_F;
final.mu_CR = mu_CR;
final.cv = cv;
final.nvc = nvc;

out = finishoutput(out, X, f, counteval, 'final', final);
end
