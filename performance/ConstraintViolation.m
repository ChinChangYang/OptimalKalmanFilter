function [c, ceq] = ConstraintViolation(x)
%CONSTRAINTVIOLATION Constraint violation measure
% x: Kalman filter (A11; A12; A22; C2)

lb = [0.3; -0.05; -0.05; 0.9];
ub = [0.5; 0.25; 0.45; 1.1];

A	= [	x(1),	x(2); ...
		-0.1,	x(3)];
C	= [ 0,		x(4)];

c = zeros(1, 2 + numel(lb));
c(1) = max(0, max(abs(eig(A))) - 1);
c(2) = length(A) - rank(obsv(A, C));

for i = 1 : numel(lb)
	c(2 + i) = max(0, lb(i) - x(i));
	c(2 + i) = max(c(2 + i), x(i) - ub(i));
end

ceq = 0;
end
