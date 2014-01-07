Max Problem:

Given a Kalman filter (A, C, Q, R):

A = [	+1.0000,	+0.0100; ...
		+0.0000,	+1.0000];

C = [	+1.0000,	+0.0000];


Q = [	+0.0001,	+0.0000; ...
		+0.0000,	+0.0001];

R =		+0.0001;

Find a system to maximize the mean square error. The results can be found 
in XXX.mat.

===========================================================================

Min-Max Problem:

Find a Kalman filter Kf and a system Sys to minimize and maximize the mean 
square error f(Kf, Sys):

	min max f(Kf, Sys)
	Kf  Sys

The results can be found in XXX.mat.

===========================================================================

Min-Min-Max Problem:

Find a weight w, a Kalman filter Kf, and a system Sys, to minimize, 
minimize, and maximize the mean square error f(w, Kf, Sys):

	min min max f(w, Kf, Sys)
	 w  Kf  Sys

The results can be found in XXX.mat.

===========================================================================

Max-Min-Max Problem:

Find a weight w, a Kalman filter Kf, and a system Sys, to maximize, 
minimize, and maximize the mean square error f(w, Kf, Sys):

	max min max f(w, Kf, Sys)
	 w  Kf  Sys

The results can be found in XXX.mat.

===========================================================================