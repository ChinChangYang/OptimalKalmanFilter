function f = error_nominal(x)
% Establish the table with Kalman Gain
% x: [Kalman filter, system]'
weight = 1;

A	= [	x(1),	x(2); ...
		0,		x(3)];
C	= [ x(4),	0];
Q	= [ 0.0001,	0; ...
		0,		0.0001];
R	= 0.0001;

K = FeedbackGain(A, C, Q, R);
f = cost_zeta(x(1:4), K, x(5:8), weight);
end
