function f = error_zeta(x, y, z)
% Establish the table with Kalman Gain
% x: weight
% y: filter
% z: system
A	= [	y(1),	y(2); ...
		0,		y(3)];
C	= [ y(4),	0];
Q	= [ 0.0001,	0; ...
		0,		0.0001];
R	= 0.0001;

K = FeedbackGain(A, C, Q, R);
f = cost_zeta(y, K, z, x);
end
