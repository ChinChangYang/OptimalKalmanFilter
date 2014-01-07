function [cost, X, X_h] = cost_zeta(K_AC, K, S_AC, weight)
% K_AC:		Kalman filter
% K:		Kalman Gain
% S_AC:		System
% weight:	Weight

A11 = S_AC(1);
A12 = S_AC(2);
A22 = S_AC(3);
C12 = S_AC(4);

%%%%% initialize
data_amount=400;
state=2;
in=2;

A_Kf = [K_AC(1)  K_AC(2); -0.1  K_AC(3)];
C_Kf = [0 K_AC(4)];
B_Kf = eye(2);

%%%%% convert the system(X) into the observe form(Xo)
A_Sys = [A11 A12; -0.1  A22];
B_Sys = eye(2);
C_Sys = [0 C12];

W = getW;
V = getV;

X			= zeros(state, data_amount);
X(:,1)		= [0.1;-0.1];
X_h			= zeros(state, data_amount);
Y			= zeros(1, data_amount);
Y(1)		= C_Sys * X(:, 1) + V(:, 1);
U			= zeros(in, data_amount);
e			= zeros(1, data_amount);
e(1)		= Y(1);
sum = (X(:, 1) - X_h(:, 1))' * (X(:, 1) - X_h(:, 1));
for k = 2 : data_amount
   X(:, k) = A_Sys * X(:, k - 1) + B_Sys * U(:, k - 1) + W(:, k - 1);
   Y(:, k) = C_Sys * X(:, k) + V(:, k);
   
   X_h(:, k) = A_Kf * X_h(:, k - 1) + B_Kf * U(:, k - 1) + K * weight * e(:,k-1);
   e(:, k) = Y(:, k) - C_Kf * X_h(:, k);
   sum = sum + (X(:, k) - X_h(:, k))' * (X(:, k) - X_h(:, k));
end

cost = sum / data_amount;