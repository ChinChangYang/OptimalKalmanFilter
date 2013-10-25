function [X_bar] = noisefree_state(S_AC)
% S_AC:		System

A11 = S_AC(1);
A12 = S_AC(2);
A22 = S_AC(3);

%%%%% initialize
data_amount=400;
state=2;
in=2;

%%%%% convert the system(X) into the observe form(Xo)
A_Sys = [A11 A12; -0.1  A22];
B_Sys = eye(2);

X_bar			= zeros(state, data_amount);
X_bar(:,1)		= [0.1;-0.1];
U			= zeros(in, data_amount);
for k = 2 : data_amount
   X_bar(:, k) = A_Sys * X_bar(:, k - 1) + B_Sys * U(:, k - 1);
end