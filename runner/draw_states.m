function draw_states
%DRAW_STATES Draw system state and its estimate

%% Min-max Problem for Kalman Filtering
load('minmax_results_by_minmaxtcjadebin.mat');

A	= [	xminmax1(1),	xminmax1(2); ...
		-0.1,			xminmax1(3)];
C	= [ 0,				xminmax1(4)];
Q	= [ 0.005,			0; ...
		0,				0.005];
R	= 0.001;

K = FeedbackGain(A, C, Q, R);
[f, X, X_h] = cost_zeta(xminmax1, K, xminmax2, 1);
X_bar = noisefree_state(xminmax2);
figure(1);
plot(X(1, :), 'b');
hold('on');
plot(X_h(1, :), 'r');
hold('on');
plot(X_bar(1, :), 'k');
hold('off');
title('Min-max Problem for Kalman Filtering');
xlabel('Time index');
ylabel('Magnitude');
fprintf('Min-max: f = %.4E\n', f);

%% Min-min-max Problem for Kalman Filtering
load('minminmax_results_by_minmaxtcjadebin.mat');

A	= [	xminmax1(1),	xminmax1(2); ...
		-0.1,			xminmax1(3)];
C	= [ 0,				xminmax1(4)];
Q	= [ 0.005,			0; ...
		0,				0.005];
R	= 0.001;

K = FeedbackGain(A, C, Q, R);
[f, X, X_h] = cost_zeta(xminmax1, K, xminmax2, xminmax1(5));
X_bar = noisefree_state(xminmax2);
figure(2);
plot(X(1, :), 'b');
hold('on');
plot(X_h(1, :), 'r');
hold('on');
plot(X_bar(1, :), 'k');
hold('off');
title('Min-min-max Problem for Kalman Filtering');
xlabel('Time index');
ylabel('Magnitude');
fprintf('Min-min-max: f = %.4E\n', f);

%% Max-min-max Problem for Kalman Filtering
load('maxminmax_results_by_maxminmaxtcjade.mat');
xminmax1 = xbest2;
xminmax2 = xbest3;

A	= [	xminmax1(1),	xminmax1(2); ...
		-0.1,			xminmax1(3)];
C	= [ 0,				xminmax1(4)];
Q	= [ 0.005,			0; ...
		0,				0.005];
R	= 0.001;

K = FeedbackGain(A, C, Q, R);
[f, X, X_h] = cost_zeta(xminmax1, K, xminmax2, xbest1);
X_bar = noisefree_state(xminmax2);
figure(3);
plot(X(1, :), 'b');
hold('on');
plot(X_h(1, :), 'r');
hold('on');
plot(X_bar(1, :), 'k');
hold('off');
title('Max-min-max Problem for Kalman Filtering');
xlabel('Time index');
ylabel('Magnitude');
fprintf('Max-min-max: f = %.4E\n', f);
end
