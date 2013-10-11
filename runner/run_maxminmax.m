function run_maxminmax
%RUN_MAXMINMAX Run the max-min-max problem
startTime = tic;
close all;
rng(1, 'twister');
solver = 'maxminmaxtcjade';
fitfun = 'error_zeta';
maxfunevals = 3e6;
solverOptions1.dimensionFactor = 5;
solverOptions1.F = 0.7;
solverOptions1.CR = 0.99;
solverOptions1.TolX = eps;
solverOptions1.TolFun = eps;
solverOptions1.TolX_DecayRate = 0.5;
solverOptions1.Display = 'iter';
solverOptions1.RecordPoint = 1000;
solverOptions2.dimensionFactor = 5;
solverOptions3.dimensionFactor = 5;
lb1 = 0.7;
ub1 = 1.3;
lb2 = [0.3; -0.05; -0.05; 0.9];
ub2 = [0.5; 0.25; 0.45; 1.1];
lb3 = lb2;
ub3 = ub2;

[xbest1, xbest2, xbest3, fbest, out] = ...
    feval(solver, fitfun, ...
	maxfunevals, lb1, ub1, lb2, ub2, lb3, ub3, solverOptions1, ...
	solverOptions2, solverOptions3);

% fprintf('FEs = \n');
% disp(out.fes(end));
fprintf('xbest1 = \n');
disp(xbest1);
fprintf('xbest2 = \n');
disp(xbest2);
fprintf('xbest3 = \n');
disp(xbest3);
fprintf('fbest = %.4E\n', fbest);
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
alternative_angle = out.angle;
alternative_angle(out.angle > pi/4) = out.angle(out.angle > pi/4) - pi/2;
if std(out.angle) < std(alternative_angle)
	semilogx(out.fes, out.angle / pi * 180, 'k');
	axis([0, out.fes(end), 0, 90]);
else
	semilogx(out.fes, alternative_angle / pi * 180, 'k');
	axis([0, out.fes(end), -45, 45]);
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

if isfield(out, 'success_rate')
	figure;
	plot(out.fes, out.success_rate);
	title(sprintf('Solve %s by %s', fitfun, solver));
	xlabel('FEs');
	ylabel('success_rate');
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

save(sprintf('maxminmax_results_by_%s.mat', solver), ...
	'xbest1', 'xbest2', 'xbest3', 'fbest', 'out');

toc(startTime);
end
