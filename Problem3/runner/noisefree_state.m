function [X_bar] = noisefree_state(S_AC)
% S_AC:		System

%%%%% initialize
data_amount=400;
state=2;

%%%%% convert the system(X) into the observe form(Xo)
A_Sys = [S_AC(1),	S_AC(2); ...
		0,			S_AC(3)];

X_bar			= zeros(state, data_amount);
X_bar(:,1)		= [1; 1];
for k = 2 : data_amount
   X_bar(:, k) = A_Sys * X_bar(:, k - 1);
end