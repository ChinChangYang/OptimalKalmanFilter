function K = FeedbackGain(A, C, Q, R)
%FEEDBACKGAIN State feedback gain
% x: Kalman filter

[~, ~, K] = dare(A', C', Q, R);
K = K';
end
