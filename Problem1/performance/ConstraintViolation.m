function [c, ceq] = ConstraintViolation(x, ~)
%CONSTRAINTVIOLATION Constraint violation measure
% x: Kalman filter (A11; A12; A22; C2)

A	= [	x(1),	x(2); ...
		-0.1,	x(3)];
C	= [ 0,		x(4)];

c = zeros(1, 2);
c(1) = max(0, max(abs(eig(A))) - 1);
c(2) = length(A) - rank(obsv(A, C));

ceq = 0;
end
