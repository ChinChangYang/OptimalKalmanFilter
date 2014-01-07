function K = FeedbackGain(A, C, Q, R)
%FEEDBACKGAIN State feedback gain
% x: Kalman filter

try
	[~, ~, K] = dare(A', C', Q, R);
catch ME1 %#ok<NASGU>
	try
		C = C + eps + max(abs(A(:))) * eps + max(abs(C(:))) * eps;
		[~, ~, K] = dare(A', C', Q, R);
	catch ME2 %#ok<NASGU>
		A = 0 * A + eps;
		C = 0 * C + eps;
		[~, ~, K] = dare(A', C', Q, R);
	end
end

K = K';
end
