function f = error_zeta(x,y,z)
% Establish the table with Kalman Gain
% x: weight
% y: filter
% z: system
[flag,K] = KalmanGain(y);
if (flag == 1)
    f = 1e10;  %penalty
    return;
end
[f,flag] = cost_zeta(y,K,z,x);
if (flag == 1)
    f = 1e10;  %penalty
end
