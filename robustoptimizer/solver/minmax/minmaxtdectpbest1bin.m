function [xminmax1, xminmax2, fminmax, out] = minmaxtdectpbest1bin(fitfun1, ...
	maxfunevals1, lb1, ub1, lb2, ub2, options1, options2)
% MINMAXTCDE Min-Max Tracer Coevolutionary Differential Evolution with
% current to pbest mutation strategy and binomial crossover operator
% MINMAXTCDE(fitfun, maxfunevals1, lb1, ub1, lb2, ub2) minimizes the
% function fitfun1 associated with a maximizer among box limitations [lb1,
% ub1] of minimizers and [lb2, ub2] of maximizers for the maximal function
% evaluations maxfunevals1.
% MINMAXTCDE(..., options1) minimizes the function with the given options
% Options1 for the 1st layer.
% MINMAXTCDE(..., options1, options2) minimize the function with the given
% options Options2 for the 2nd layer.
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
defaultOptions1.mu_CR = 0.5;
defaultOptions1.mu_F = 1.0;
defaultOptions1.delta_CR = 0.1;
defaultOptions1.delta_F = 0.1;
defaultOptions1.p = 0.05;
defaultOptions1.w = 0.1;
defaultOptions1.Display = 'off';
defaultOptions1.FactorNP = 1;
defaultOptions1.Restart = 0;
defaultOptions1.RecordPoint = 100;
defaultOptions1.Noise = false;
defaultOptions1.SampleFactor = 1.01;
defaultOptions1.ftarget = -Inf;
defaultOptions1.TolFun = eps;
defaultOptions1.TolX = 100 * eps;
options1 = setdefoptions(options1, defaultOptions1);

% Default options for Layer 2
defaultOptions2 = defaultOptions1;
defaultOptions2.Display = 'off';
defaultOptions2.RecordPoint = 1;
defaultOptions2.TolX = 0;
defaultOptions2.TolFun = 0;
options2 = setdefoptions(options2, defaultOptions2);
maxfunevals2 = 2 * ceil(options2.dimensionFactor * D2);

% Initialize algorithmic variables
dimensionFactor = max(1, options1.dimensionFactor);
mu_CR = options1.mu_CR;
mu_F = options1.mu_F;
delta_CR = options1.delta_CR;
delta_F = options1.delta_F;
p = options1.p;
w = options1.w;
isDisplayIter = strcmp(options1.Display, 'iter');
RecordPoint = max(1, floor(options1.RecordPoint));
TolFun = options1.TolFun;
TolX = options1.TolX;

NP1 = ceil(dimensionFactor * D1);
NP2 = ceil(options2.dimensionFactor * D2);
FEsPerIter1 = NP1 * maxfunevals2;

% Initialize variables
counteval = 0;
countiter = 1;
out = initoutput(RecordPoint, D1, NP1, maxfunevals1, ...
	'mu_F', 'mu_CR');

% Initialize contour data
if isDisplayIter
	[XX, YY, ZZ] = minmaxcontourdata(D1, lb1, ub1, lb2, ub2, fitfun1);
end

% Initialize population
if NP1 < 1e1
	LHS = lhsdesign(NP1, D1, 'iteration', 10)';
elseif NP1 < 1e2
	LHS = lhsdesign(NP1, D1, 'iteration', 2)';
else
	LHS = rand(D1, NP1);
end

X1 = zeros(D1, NP1);
for i = 1 : NP1
	X1(:, i) = lb1 + (ub1 - lb1) .* LHS(:, i);
end

% Initialize variables
Xmax2 = zeros(D2, NP1);
Umax2 = Xmax2;
V1 = X1;
U1 = X1;
X2 = zeros(D2, NP2, NP1);
f2 = zeros(NP2, NP1);
mu_F2 = zeros(1, NP1);
mu_CR2 = zeros(1, NP1);
V2 = X2;
pbest_size = p * NP1;

% Evaluation
f = zeros(1, NP1);
for i = 1 : NP1
	fitfunX2i = @(y) -feval(fitfun1, X1(:, i), y);
	[Xmax2(:, i), fmax, out2] = jadebin(fitfunX2i, lb2, ub2, maxfunevals2, options2);
	X2(:, :, i) = out2.final.X;
	f2(:, i) = out2.final.f;
	mu_F2(i) = out2.final.mu_F;
	mu_CR2(i) = out2.final.mu_CR;
	counteval = counteval + out2.fes(end);
	f(i) = -fmax;
end

% Sort
[f, fidx] = sort(f);
X1 = X1(:, fidx);
X2 = X2(:, :, fidx);
f2 = f2(:, fidx);
mu_F2 = mu_F2(fidx);
mu_CR2 = mu_CR2(fidx);

% Display
if isDisplayIter
	displayitermessages([X1; Xmax2], [X1; Xmax2], f, countiter, XX, YY, ZZ, 'mu_F', mu_F, 'mu_CR', mu_CR);
end

% Record minimal function values
out = updateoutput(out, X1, f, counteval, ...
	'mu_F', mu_F, 'mu_CR', mu_CR);
countiter = countiter + 1;

while true
	% Termination conditions	
	stdf = std(f);
	stdX1 = std(X1, 0, 2);
	meanX1 = mean(X1, 2);
	outofmaxfunevals = counteval > maxfunevals1 - FEsPerIter1;
	fitnessconvergence = stdf <= mean(f) * 100 * eps || stdf <= realmin || stdf <= TolFun;
	solutionconvergence = all(stdX1 <= meanX1 * 100 * eps) || all(stdX1 <= 100 * realmin) || ...
			all(stdX1 <= TolX);
	
	% Convergence conditions	
	if outofmaxfunevals || fitnessconvergence || solutionconvergence
		break;
	end
	
	% Scaling factor and crossover rate
	S_F = zeros(1, NP1);
	S_CR = zeros(1, NP1);
	CR = mu_CR + delta_CR * randn(1, NP1);
	CR(CR > 1) = 1;
	CR(CR < 0) = 0;
	F = cauchyrnd(mu_F, delta_F, NP1, 1);
	F(F > 1) = 1 - eps - 0.01 * rand;
	
	for retry = 1 : 3
		if all(F > 0)
			break;
		end
		
		F(F <= 0) = cauchyrnd(mu_F, delta_F, sum(F <= 0), 1);
		F(F > 1) = 1 - eps - 0.01 * rand;
	end
	
	F(F <= 0) = 100 * eps * (1 + rand);
	
	Succ_Counter = 0;
	
	% Mutation
	for i = 1 : NP1
		
		for checkRetry = 1 : 3
			
			% Generate pbest_idx
			for retry = 1 : 3
				pbest_idx = max(1, ceil(rand * pbest_size));
				if ~all(X1(:, pbest_idx) == X1(:, i))
					break;
				end
			end
			
			% Generate r1
			for retry = 1 : 3
				r1 = floor(1 + NP1 * rand);
				if i ~= r1
					break;
				end
			end
			
			% Generate r2
			for retry = 1 : 3
				r2 = floor(1 + NP1 * rand);
				if ~(all(X1(:, i) == X1(:, r2)) || all(X1(:, r1) == X1(:, r2)))
					break;
				end
			end
			
			V1(:, i) = X1(:, i) + F(i) * (X1(:, pbest_idx) - X1(:, i) + ...
				X1(:, r1) - X1(:, r2));
			
			V2(:, :, i) = X2(:, :, i) + F(i) * (X2(:, :, pbest_idx) - X2(:, :, i) + ...
				X2(:, :, r1) - X2(:, :, r2));
			
			% Check boundary
			if all(V1(:, i) > lb1) && all(V1(:, i) < ub1)
				break;
			end
		end
	end
	
	for i = 1 : NP1
		% Binominal Crossover
		jrand = floor(1 + D1 * rand);
		for j = 1 : D1
			if rand < CR(i) || j == jrand
				U1(j, i) = V1(j, i);
			else
				U1(j, i) = X1(j, i);
			end
		end
	end
	
	% Repair
	for i = 1 : NP1
		for j = 1 : D1
			if U1(j, i) < lb1(j)
				U1(j, i) = 2 * lb1(j) - U1(j, i);
			elseif U1(j, i) > ub1(j)
				U1(j, i) = 2 * ub1(j) - U1(j, i);
			else
				continue;
			end
			
			if U1(j, i) < lb1(j)
				U1(j, i) = lb1(j);
			elseif U1(j, i) > ub1(j)
				U1(j, i) = ub1(j);
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
	
	% Selection
	for i = 1 : NP1
		% Compute fxi, f(i)
		fitfunX2i = @(X2) -feval(fitfun1, X1(:, i), X2);
		optionsX2i = options2;
		optionsX2i.initial.X = X2(:, :, i);
		optionsX2i.initial.f = f2(:, i);
		optionsX2i.initial.mu_F = mu_F2(i);
		optionsX2i.initial.mu_CR = mu_CR2(i);
		[Xmax2(:, i), fmax, outX2i] = jadebin(fitfunX2i, lb2, ub2, maxfunevals2, optionsX2i);
		counteval = counteval + outX2i.fes(end);
		f(i) = -fmax;
		
		% Compute fui
		fitfunU2i = @(U2) -feval(fitfun1, U1(:, i), U2);
		optionsU2i = options2;
		optionsU2i.initial.X = V2(:, :, i);
		optionsU2i.initial.f = [];
		optionsU2i.initial.mu_F = mu_F2(i);
		optionsU2i.initial.mu_CR = mu_CR2(i);
		[Umax2(:, i), fmax, outU2i] = jadebin(fitfunU2i, lb2, ub2, countiter * maxfunevals2, optionsU2i);
		counteval = counteval + outU2i.fes(end);
		fui = -fmax;
		
		% Replacement
		if fui < f(i)
			f(i) = fui;
			X1(:, i) = U1(:, i);
			Xmax2(:, i) = Umax2(:, i);
			S_CR(Succ_Counter + 1) = CR(i);
			S_F(Succ_Counter + 1) = F(i);
			Succ_Counter = Succ_Counter + 1;			
			
			X2(:, :, i) = outU2i.final.X;
			f2(:, i) = outU2i.final.f;
			mu_F2(i) = outU2i.final.mu_F;
			mu_CR2(i) = outU2i.final.mu_CR;
		else			
			X2(:, :, i) = outX2i.final.X;
			f2(:, i) = outX2i.final.f;
			mu_F2(i) = outX2i.final.mu_F;
			mu_CR2(i) = outX2i.final.mu_CR;
		end
	end
	
	% Display
	if isDisplayIter
		displayitermessages(...
			[X1; Xmax2], [U1; Umax2], f, countiter, XX, YY, ZZ, ...
			'mu_F', mu_F, 'mu_CR', mu_CR);
	end
	
	% Update CR and F
	if Succ_Counter > 0
		mu_CR = (1-w) * mu_CR + w * mean(S_CR(1 : Succ_Counter));
		mu_F = (1-w) * mu_F + w * sum(S_F(1 : Succ_Counter).^2) / sum(S_F(1 : Succ_Counter));
	end
	
	% Sort
	[f, fidx] = sort(f);
	X1 = X1(:, fidx);
	X2 = X2(:, :, fidx);
	f2 = f2(:, fidx);
	mu_F2 = mu_F2(fidx);
	mu_CR2 = mu_CR2(fidx);
	
	% Record
	out = updateoutput(out, X1, f, counteval, ...
		'mu_F', mu_F, 'mu_CR', mu_CR);
	countiter = countiter + 1;
end

[fminmax, fminidx] = min(f);
xminmax1 = X1(:, fminidx);
xminmax2 = Xmax2(:, fminidx);
out = finishoutput(out, X1, f, counteval);
end