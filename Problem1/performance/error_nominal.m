function f = error_nominal(x)
% Establish the table with Kalman Gain
% x: [Kalman filter, system]'
weight = 1;

A	= [	x(1),	x(2); ...
		-0.1,	x(3)];
C	= [ 0,		x(4)];
Q	= [ 0.005,	0; ...
		0,		0.005];
R	= 0.001;

K = FeedbackGain(A, C, Q, R);
f = cost_zeta(x(1:4), K, x(5:8), weight);
end
