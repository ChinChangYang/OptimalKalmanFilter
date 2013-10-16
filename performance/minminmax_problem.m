function f = minminmax_problem(x, y)
% Establish the table with Kalman Gain
% x: [filter; weight]
% y: system
A	= [	x(1),	x(2); ...
		-0.1,	x(3)];
	
C	= [ 0,		x(4)];

Q	= [ 0.005,	0; ...
		0,		0.005];
	
R	= 0.001;

K = FeedbackGain(A, C, Q, R);
f = cost_zeta(x, K, y, x(5));
end
