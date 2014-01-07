function f = error_nominal_zeta(x)
% Establish the table with Kalman Gain
% x: [Kalman filter, system, weight]'
A	= [	x(1),	x(2); ...
		-0.1,	x(3)];
C	= [ 0,		x(4)];
Q	= [ 0.005,	0; ...
		0,		0.005];
R	= 0.001;

K = FeedbackGain(A, C, Q, R);
f = cost_zeta(x(1:4), K, x(5:8), x(9));
end
