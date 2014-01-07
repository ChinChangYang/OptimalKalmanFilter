function estimate_time
%ESTIMATE_TIME Estimate the running of Kalman filtering
startTime = tic;
T = 100;
lb1 = 0.7;
ub1 = 1.3;
lb2 = [0.3; -0.05; -0.05; 0.9];
ub2 = [0.5; 0.25; 0.45; 1.1];
lb3 = lb2;
ub3 = ub2;
D1 = numel(lb1);
D2 = numel(lb2);
D3 = numel(lb3);
f = zeros(1, T);

for t = 1 : T
	x = lb1 + rand(D1, 1) .* (ub1 - lb1);
	y = lb2 + rand(D2, 1) .* (ub2 - lb2);
	z = lb3 + rand(D3, 1) .* (ub3 - lb3);
	f(t) = error_zeta(x, y, z);
end

toc(startTime);
end

