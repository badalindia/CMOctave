function [Q Qx] = GALS1d(x, Q, Qx, U, t0, tf, dt, ifPlot)

% the function perform time stepping of advection equation using 'Gradient Augmented Level Sets' method. 
% x is the spatial domain
% Q is initial field on x
% Qx is gradient of Q
% U is velocity  (function handle @(X, T))
% t0: initial time
% tf: final time
% ifPlot: (boolean) if data has to be plotted for each time step

h = x(2) - x(1);
t = t0;

I = length(x);
xMin = x(1);
xMax = x(I);
xDelta = xMax - xMin;

while t<tf

	if(t+dt>tf)
		dt = tf - t;
	end;
	
	%trick to find derivative
	ex = 10^(-6);
	x0p = RK3(x+ex, t+dt, -dt, U);
	x0n = RK3(x-ex, t+dt, -dt, U);
	
	%adding warping
	x0p = x0p + ((x0p<xMin) - (x0p>=xMax))*xDelta;
	x0n = x0n + ((x0n<xMin) - (x0n>=xMax))*xDelta;

	Q0p = H1d(x, Q, Qx, x0p);
	Q0n = H1d(x, Q, Qx, x0n);

	Q  = (Q0p + Q0n)/(2);
	Qx = (Q0p - Q0n)/(2*ex);

	t = t+dt;
	
	if(ifPlot)
		plot(x, Q, '-or');
		axis([xMin, xMax, -1.5, 1.5]);
		drawnow;
	end

end


