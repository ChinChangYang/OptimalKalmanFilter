function draw_states
%DRAW_STATES Draw system state and its estimate
close all;
%% Min-min Problem for Kalman Filtering
load('minmin_results_by_jadebin_201311081107.mat');

A	= [	xmin(1),	xmin(2); ...
		0,			xmin(3)];
C	= [ xmin(4),	0];
Q	= [ 0.0001,		0; ...
		0,			0.0001];
R	= 0.0001;

K = FeedbackGain(A, C, Q, R);
[f, X, X_h] = cost_zeta(xmin(1:4), K, xmin(5:8), 1);
X_bar = noisefree_state(xmin(5:8));
figure(1);
plot(X(1, :), 'b');
hold('on');
plot(X_h(1, :), 'r');
hold('on');
plot(X_bar(1, :), 'k');
hold('off');
title('Min-min Problem for Kalman Filtering');
xlabel('Time index');
ylabel('Magnitude');
print('state_01.tif', '-dtiff', '-r400');
figure(2);
plot(X(2, :), 'b');
hold('on');
plot(X_h(2, :), 'r');
hold('on');
plot(X_bar(2, :), 'k');
hold('off');
title('Min-min Problem for Kalman Filtering');
xlabel('Time index');
ylabel('Magnitude');
print('state_02.tif', '-dtiff', '-r400');
fprintf('Min-min: f = %.4E\n', f);

%% Max-max Problem for Kalman Filtering
load('maxmax_results_by_jadebin_201311081111.mat');

A	= [	xmin(1),	xmin(2); ...
		0,			xmin(3)];
C	= [ xmin(4),	0];
Q	= [ 0.0001,		0; ...
		0,			0.0001];
R	= 0.0001;

K = FeedbackGain(A, C, Q, R);
[f, X, X_h] = cost_zeta(xmin(1:4), K, xmin(5:8), 1);
X_bar = noisefree_state(xmin(5:8));
figure(3);
plot(X(1, :), 'b');
hold('on');
plot(X_h(1, :), 'r');
hold('on');
plot(X_bar(1, :), 'k');
hold('off');
title('Max-max Problem for Kalman Filtering');
xlabel('Time index');
ylabel('Magnitude');
print('state_03.tif', '-dtiff', '-r400');
figure(4);
plot(X(2, :), 'b');
hold('on');
plot(X_h(2, :), 'r');
hold('on');
plot(X_bar(2, :), 'k');
hold('off');
title('Max-max Problem for Kalman Filtering');
xlabel('Time index');
ylabel('Magnitude');
print('state_04.tif', '-dtiff', '-r400');
fprintf('Max-max: f = %.4E\n', f);

%% Min-max Problem for Kalman Filtering
load('minmax_results_by_mmdeb1b_pce_201401021315.mat');

A	= [	xminmax1(1),	xminmax1(2); ...
		0,				xminmax1(3)];
C	= [ xminmax1(4),	0];
Q	= [ 0.0001,			0; ...
		0,				0.0001];
R	= 0.0001;

K = FeedbackGain(A, C, Q, R);
[f, X, X_h] = cost_zeta(xminmax1, K, xminmax2, 1);
X_bar = noisefree_state(xminmax2);
figure(5);
plot(X(1, :), 'b');
hold('on');
plot(X_h(1, :), 'r');
hold('on');
plot(X_bar(1, :), 'k');
hold('off');
title('Min-max Problem for Kalman Filtering');
xlabel('Time index');
ylabel('Magnitude');
print('state_05.tif', '-dtiff', '-r400');
figure(6);
plot(X(2, :), 'b');
hold('on');
plot(X_h(2, :), 'r');
hold('on');
plot(X_bar(2, :), 'k');
hold('off');
title('Min-max Problem for Kalman Filtering');
xlabel('Time index');
ylabel('Magnitude');
print('state_06.tif', '-dtiff', '-r400');
fprintf('Min-max: f = %.4E\n', f);

%% Min-min-min Problem for Kalman Filtering
load('minminmin_results_by_jadebin_201311081236.mat');

A	= [	xmin(1),	xmin(2); ...
		0,			xmin(3)];
C	= [ xmin(4),	0];
Q	= [ 0.0001,		0; ...
		0,			0.0001];
R	= 0.0001;

K = FeedbackGain(A, C, Q, R);
[f, X, X_h] = cost_zeta(xmin(1:4), K, xmin(5:8), 1);
X_bar = noisefree_state(xmin(5:8));
figure(7);
plot(X(1, :), 'b');
hold('on');
plot(X_h(1, :), 'r');
hold('on');
plot(X_bar(1, :), 'k');
hold('off');
title('Min-min-min Problem for Kalman Filtering');
xlabel('Time index');
ylabel('Magnitude');
print('state_07.tif', '-dtiff', '-r400');
figure(8);
plot(X(2, :), 'b');
hold('on');
plot(X_h(2, :), 'r');
hold('on');
plot(X_bar(2, :), 'k');
hold('off');
title('Min-min-min Problem for Kalman Filtering');
xlabel('Time index');
ylabel('Magnitude');
print('state_08.tif', '-dtiff', '-r400');
fprintf('Min-min-min: f = %.4E\n', f);

%% Max-max-max Problem for Kalman Filtering
load('maxminmax_results_by_mmmdeb1b_pce_201401021930.mat');

A	= [	xmin(1),	xmin(2); ...
		0,			xmin(3)];
C	= [ xmin(4),	0];
Q	= [ 0.0001,		0; ...
		0,			0.0001];
R	= 0.0001;

K = FeedbackGain(A, C, Q, R);
[f, X, X_h] = cost_zeta(xmin(1:4), K, xmin(5:8), 1);
X_bar = noisefree_state(xmin(5:8));
figure(9);
plot(X(1, :), 'b');
hold('on');
plot(X_h(1, :), 'r');
hold('on');
plot(X_bar(1, :), 'k');
hold('off');
title('Max-max-max Problem for Kalman Filtering');
xlabel('Time index');
ylabel('Magnitude');
print('state_09.tif', '-dtiff', '-r400');
figure(10);
plot(X(2, :), 'b');
hold('on');
plot(X_h(2, :), 'r');
hold('on');
plot(X_bar(2, :), 'k');
hold('off');
title('Max-max-max Problem for Kalman Filtering');
xlabel('Time index');
ylabel('Magnitude');
print('state_10.tif', '-dtiff', '-r400');
fprintf('Max-max-max: f = %.4E\n', f);

%% Min-min-max Problem for Kalman Filtering
load('minminmax_results_by_mmdeb1b_pce_201401021352.mat');

A	= [	xminmax1(1),	xminmax1(2); ...
		0,				xminmax1(3)];
C	= [ xminmax1(4),	0];
Q	= [ 0.0001,			0; ...
		0,				0.0001];
R	= 0.0001;

K = FeedbackGain(A, C, Q, R);
[f, X, X_h] = cost_zeta(xminmax1, K, xminmax2, xminmax1(5));
X_bar = noisefree_state(xminmax2);
figure(11);
plot(X(1, :), 'b');
hold('on');
plot(X_h(1, :), 'r');
hold('on');
plot(X_bar(1, :), 'k');
hold('off');
title('Min-min-max Problem for Kalman Filtering');
xlabel('Time index');
ylabel('Magnitude');
print('state_11.tif', '-dtiff', '-r400');
figure(12);
plot(X(2, :), 'b');
hold('on');
plot(X_h(2, :), 'r');
hold('on');
plot(X_bar(2, :), 'k');
hold('off');
title('Min-min-max Problem for Kalman Filtering');
xlabel('Time index');
ylabel('Magnitude');
print('state_12.tif', '-dtiff', '-r400');
fprintf('Min-min-max: f = %.4E\n', f);

%% Max-min-max Problem for Kalman Filtering
load('maxminmax_results_by_mmmdeb1b_pce_201401021930.mat');

A	= [	xbest2(1),		xbest2(2); ...
		0,				xbest2(3)];
C	= [ xbest2(4),		0];
Q	= [ 0.0001,			0; ...
		0,				0.0001];
R	= 0.0001;

K = FeedbackGain(A, C, Q, R);
[f, X, X_h] = cost_zeta(xbest2, K, xbest3, xbest1);
X_bar = noisefree_state(xbest3);
figure(13);
plot(X(1, :), 'b');
hold('on');
plot(X_h(1, :), 'r');
hold('on');
plot(X_bar(1, :), 'k');
hold('off');
title('Max-min-max Problem for Kalman Filtering');
xlabel('Time index');
ylabel('Magnitude');
print('state_13.tif', '-dtiff', '-r400');
figure(14);
plot(X(2, :), 'b');
hold('on');
plot(X_h(2, :), 'r');
hold('on');
plot(X_bar(2, :), 'k');
hold('off');
title('Max-min-max Problem for Kalman Filtering');
xlabel('Time index');
ylabel('Magnitude');
print('state_14.tif', '-dtiff', '-r400');
fprintf('Max-min-max: f = %.4E\n', f);
end
