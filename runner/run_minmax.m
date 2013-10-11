function [xminmax1, xminmax2, fminmax, out] = run_minmax
%RUN_MINMAX Run the minimax problem
startTime = tic;
close all;
rng(1, 'twister');
solver = 'minmaxtcdebest1bin';
fitfun = 'error_nozeta';
D_Min = 4;
D_Max = 4;
D = D_Min + D_Max;
maxfunevals = 1e5;
solverOptions1.dimensionFactor = 5;
solverOptions1.F = 0.7;
solverOptions1.CR = 0.99;
solverOptions1.delta_CR = 0.01;
solverOptions1.TolX = 1e-11;
solverOptions1.TolFun = 0;
solverOptions1.Display = 'iter';
solverOptions1.RecordPoint = 1000;
solverOptions1.TolX_DecayRate = 0.5;
solverOptions2.dimensionFactor = 5;
solverOptions2.F = 0.7;
solverOptions2.CR = 0.5;
solverOptions2.TolX = 0;
solverOptions2.TolFun = 0;
lb1 = [0.3; -0.05; -0.05; 0.9];
ub1 = [0.5; 0.25; 0.45; 1.1];
lb2 = lb1;
ub2 = ub1;
[xminmax1, xminmax2, fminmax, out] = ...
    feval(solver, fitfun, ...
	maxfunevals, lb1, ub1, lb2, ub2, solverOptions1, solverOptions2);
% fprintf('FEs = \n');
% disp(out.fes(end));
fprintf('xminmax1 = \n');
disp(xminmax1);
fprintf('xminmax2 = \n');
disp(xminmax2);
fprintf('fminmax = %.4E\n', fminmax);
figure(2);
subplot(241);
hold off;
plot(out.fes, out.fmean);
xlabel('FEs');
title('function values');
subplot(242);
semilogy(out.fes, out.fstd);
xlabel('FEs');
title('Std. of function values');
subplot(243);
hold off;
for i = 1 : numel(out.xmean(:, 1))
    plot(out.fes, out.xmean(i, :), getlinespec(i));
    hold on;
end
xlabel('FEs');
title('Mean of X solutions');
subplot(245);
hold off;
for i = 1 : numel(out.xstd(:, 1))
	semilogy(out.fes, out.xstd(i, :), getlinespec(i));
    hold on;
end
xlabel('FEs');
title('Std. of X solutions');
subplot(246);
hold off;
semilogy(out.fes, out.distancemin, 'r');
hold on;
semilogy(out.fes, out.distancemax, 'r');
semilogy(out.fes, out.distancemean, 'b');
semilogy(out.fes, out.distancemedian, 'g');
xlabel('FEs');
title('Distances between each pair of X');
subplot(247);
hold off;
semilogy(out.fes, out.cond);
xlabel('FEs');
title('Condition number');
subplot(244);
hold off;
semilogy(out.fes, out.innerFstd);
xlabel('FEs');
title('Std. f in inner states');
subplot(248);
hold off;
semilogy(out.fes, out.innerMeanXstd);
xlabel('FEs');
title('Mean of std. X in inner states');

figure;
NP = D * solverOptions1.dimensionFactor;
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

if isfield(out, 'alpha')
	figure;
	plot(out.fes, out.alpha);
	title(sprintf('Solve %s by %s', fitfun, solver));
	xlabel('FEs');
	ylabel('alpha');
end

if isfield(out, 'successRate')
	figure;
	plot(out.fes, out.successRate);
	title(sprintf('Solve %s by %s', fitfun, solver));
	xlabel('FEs');
	ylabel('Success Rate');
end

if isfield(out, 'X_Converged_FEs')
	figure;
	plot(out.fes, out.X_Converged_FEs);
	title(sprintf('Solve %s by %s', fitfun, solver));
	xlabel('FEs');
	ylabel('X_Converged_FEs');
end

if isfield(out, 'U_Converged_FEs')
	figure;
	plot(out.fes, out.U_Converged_FEs);
	title(sprintf('Solve %s by %s', fitfun, solver));
	xlabel('FEs');
	ylabel('U_Converged_FEs');
end

save(sprintf('minmax_results_by_%s.mat', solver), ...
	'xminmax1', 'xminmax2', 'fminmax');

toc(startTime);
end
