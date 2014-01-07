function f = error_zeta(x, y, z)
% Establish the table with Kalman Gain
% x: weight
% y: filter
% z: system
A	= [	y(1),	y(2); ...
		-0.1,	y(3)];
C	= [ 0,		y(4)];
Q	= [ 0.005,	0; ...
		0,		0.005];
R	= 0.001;

K = FeedbackGain(A, C, Q, R);
f = cost_zeta(y, K, z, x);
end
