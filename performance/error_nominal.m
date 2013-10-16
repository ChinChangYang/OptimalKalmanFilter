function f = error_nominal(x)
% Establish the table with Kalman Gain
% x: system
nominalFilter = [0.4; 0.1; 0.2; 1.0];
weight = 1;
K = [0.0689016331734661;0.170855429044690;];
f = cost_zeta(nominalFilter, K, x, weight);
end
