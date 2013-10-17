Given a Kalman filter (A, C, Q, R):

A = [	[0.1, 0.9],		[0.1, 0.9]; ...
		-0.1000,		[0.1, 0.9]];

C = [	+0.0000,		[0.1, 0.9]];


Q = [	+0.0050,	+0.0000; ...
		+0.0000,	+0.0050];

R =		+0.0010;

===========================================================================

Min-Max Problem:

Find a Kalman filter Kf and a system Sys to minimize and maximize the mean 
square error f(Kf, Sys):

	min max f(Kf, Sys)
	Kf  Sys

The results can be found in minmax_results_by_minmaxtcjadebin_20131016.mat.

===========================================================================

Min-Min-Max Problem:

Find a weight w, a Kalman filter Kf, and a system Sys, to minimize, 
minimize, and maximize the mean square error f(w, Kf, Sys):

	min min max f(w, Kf, Sys)
	 w  Kf  Sys

The results can be found in 
minminmax_results_by_minmaxtcjadebin_20131016.mat.

Note: w in [0.7, 1.3]

===========================================================================

Max-Min-Max Problem:

Find a weight w, a Kalman filter Kf, and a system Sys, to maximize, 
minimize, and maximize the mean square error f(w, Kf, Sys):

	max min max f(w, Kf, Sys)
	 w  Kf  Sys

The results can be found in 
maxminmax_results_by_maxminmaxtcjade_20131016.mat.

Note: w in [0.7, 1.3]

===========================================================================