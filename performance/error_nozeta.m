function f = error_nozeta(x,y)
% Establish the table with Kalman Gain
% x: filter
% y: system
pert = 1;
[flag,K] = KalmanGain(x);
if (flag == 1)
    f = 1e10;  %penalty
    return;
end
[f,flag] = cost_zeta(x,K,y,pert);
if (flag == 1)
    f = 1e10;  %penalty
end