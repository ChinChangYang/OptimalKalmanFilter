function f = minminmax_problem(x, y)
% Establish the table with Kalman Gain
% x: [filter; weight]
% y: system
A	= [	x(1),	x(2); ...
		0,		x(3)];
C	= [ x(4),	0];
Q	= [ 0.0001,	0; ...
		0,		0.0001];
R	= 0.0001;

K = FeedbackGain(A, C, Q, R);
f = cost_zeta(x, K, y, x(5));
end
