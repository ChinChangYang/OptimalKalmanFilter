function run_maxmax
%RUN_MAX Run the max-max problem
startTime = tic;
close all;
rng(1, 'twister');
solver = 'jadebin';
fitfun = 'error_nominal';
fhandle = @(x) -error_nominal(x);
D = 8;
maxfunevals = 3e6;
solverOptions.dimensionFactor = 20;
solverOptions.F = 0.9;
solverOptions.CR = 0.5;
solverOptions.TolX = 0;
solverOptions.TolFun = 0;
solverOptions.Display = 'iter';
solverOptions.RecordPoint = 1000;
% solverOptions.nonlcon = 'ConstraintViolation';
lb = [0.3; -0.05; -0.05; 0.9; 0.3; -0.05; -0.05; 0.9];
ub = [0.5; 0.25; 0.45; 1.1; 0.5; 0.25; 0.45; 1.1];
[xmin, fmin, out] = ...
	feval(solver, fhandle, lb, ub, maxfunevals, solverOptions);
fprintf('out.bestever.xmin = \n');
disp(out.bestever.xmin);
if isfield(out, 'xmean')
	fprintf('xmean = \n');
	disp(out.xmean(:, end));
end
fprintf('out.bestever.fmin = \n');
disp(-out.bestever.fmin);
fprintf('xmin = \n');
disp(xmin);
fprintf('fmin = %.4E\n', -fmin);
if strncmp('bbob12', fitfun, 6)
	fprintf('fmin - fopt = %.4E\n', fmin - feval(fitfun, 'xopt'));
end
figure(2);
subplot(231);
hold off;
lowerBoundF = min([out.fmean, out.fmin]);
semilogy(out.fes, out.fmean - lowerBoundF, 'b');
hold on;
semilogy(out.fes, out.fmin - lowerBoundF, 'r');
xlabel('FEs');
title('function values');
subplot(232);
semilogy(out.fes, out.fstd);
xlabel('FEs');
title('Std. of function values');
subplot(233);
hold off;
for i = 1 : numel(out.xmean(:, 1))
	plot(out.fes, out.xmean(i, :), getlinespec(i));
	hold on;
end
xlabel('FEs');
title('Mean of X solutions');
subplot(234);
hold off;
for i = 1 : numel(out.xstd(:, 1))
	semilogy(out.fes, out.xstd(i, :), getlinespec(i));
	hold on;
end
xlabel('FEs');
title('Std. of X solutions');
subplot(235);
hold off;
semilogy(out.fes, out.distancemin, 'r');
hold on;
semilogy(out.fes, out.distancemax, 'r');
semilogy(out.fes, out.distancemean, 'b');
semilogy(out.fes, out.distancemedian, 'g');
xlabel('FEs');
title('Distances between each pair of X');
subplot(236);
hold off;
semilogy(out.fes, out.cond);
xlabel('FEs');
title('Condition number');

figure;
NP = D * solverOptions.dimensionFactor;
Gmax = maxfunevals / NP;
alternative_angle = out.angle;
alternative_angle(out.angle > pi/4) = out.angle(out.angle > pi/4) - pi/2;
if std(out.angle) < std(alternative_angle)
	semilogx(out.fes / NP, out.angle / pi * 180, 'k');
	axis([0, Gmax, 0, 90]);
else
	semilogx(out.fes / NP, alternative_angle / pi * 180, 'k');
	axis([0, Gmax, -45, 45]);
end
xlabel('Generation');
ylabel('Angle (degree)');
title('Angle between the natural basis and the eigenvector basis');

if isfield(out, 'mu_F')
	figure;
	semilogy(out.fes, out.mu_F);
	title(sprintf('Solve %s by %s', fitfun, solver));
	xlabel('FEs');
	ylabel('mu_F');
end

if isfield(out, 'mu_CR')
	figure;
	plot(out.fes, out.mu_CR);
	title(sprintf('Solve %s by %s', fitfun, solver));
	xlabel('FEs');
	ylabel('mu_CR');
end

if isfield(out, 'mu_R')
	figure;
	plot(out.fes, out.mu_R);
	title(sprintf('Solve %s by %s', fitfun, solver));
	xlabel('FEs');
	ylabel('mu_R');
end

if isfield(out, 'mu_G')
	figure;
	plot(out.fes, out.mu_G);
	title(sprintf('Solve %s by %s', fitfun, solver));
	xlabel('FEs');
	ylabel('mu_G');
end

save(sprintf('max_results_by_%s.mat', solver), ...
	'xmin', 'fmin');

toc(startTime);
end
