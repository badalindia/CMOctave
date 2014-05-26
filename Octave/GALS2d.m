function Q = GALS2d(x, y, Q, Qx, Qy, Qxy, U, V, t0, tf, dt, ifPlot)

% the function perform time stepping of advection equation using 'Gradient Augmented Level Sets' method. 
% x & y is the spatial domain (grid)
% Q, Qx, Qy, Qxy is initial field on grid
% U & V: velocity  (function handle @(X, Y, T))
% t0: initial time
% tf: final time
% ifPlot: (boolean) if data has to be plotted for each time step

page_output_immediately(true);

h = x(1, 2) - x(1, 1);
t = t0;

[Iy Ix] = size(x)
xMin = x(1, 1);
xMax = x(1, Ix);
xDelta = xMax - xMin;
yMin = y(1, 1);
yMax = y(Iy, 1);
yDelta = yMax - yMin;

ctr = 0;
if(ifPlot)
	contourf(x, y, Q, 10);
	axis([xMin, xMax, yMin, yMax]);
	drawnow;
	print( strcat('Images/Image_',num2str(Ix), '_', num2str(ctr),'.png'));
	ctr = ctr+1;
end;


while t<tf;

	if(t+dt>tf)
		dt = tf - t;
	end;

	%backward integration
	eh = 10^-6;
	[x0_xpyp, y0_xpyp] = RK3_2D(x+eh, y+eh, t+dt, -dt, U, V);
	[x0_xpyn, y0_xpyn] = RK3_2D(x+eh, y-eh, t+dt, -dt, U, V);
	[x0_xnyp, y0_xnyp] = RK3_2D(x-eh, y+eh, t+dt, -dt, U, V);
	[x0_xnyn, y0_xnyn] = RK3_2D(x-eh, y-eh, t+dt, -dt, U, V);
	
	x0_xpyp  = x0_xpyp  + ((x0_xpyp <xMin) - (x0_xpyp >=xMax))*xDelta;
	x0_xpyn  = x0_xpyn  + ((x0_xpyn <xMin) - (x0_xpyn >=xMax))*xDelta;
	x0_xnyp  = x0_xnyp  + ((x0_xnyp <xMin) - (x0_xnyp >=xMax))*xDelta;
	x0_xnyn  = x0_xnyn  + ((x0_xnyn <xMin) - (x0_xnyn >=xMax))*xDelta;
	y0_xpyp  = y0_xpyp  + ((y0_xpyp <yMin) - (y0_xpyp >=yMax))*yDelta;
	y0_xpyn  = y0_xpyn  + ((y0_xpyn <yMin) - (y0_xpyn >=yMax))*yDelta;
	y0_xnyp  = y0_xnyp  + ((y0_xnyp <yMin) - (y0_xnyp >=yMax))*yDelta;
	y0_xnyn  = y0_xnyn  + ((y0_xnyn <yMin) - (y0_xnyn >=yMax))*yDelta;
	
	Q0_xpyp  = H2dw(x, y, Q, Qx, Qy, Qxy, x0_xpyp, y0_xpyp);
	Q0_xpyn  = H2dw(x, y, Q, Qx, Qy, Qxy, x0_xpyn, y0_xpyn);
	Q0_xnyp  = H2dw(x, y, Q, Qx, Qy, Qxy, x0_xnyp, y0_xnyp);
	Q0_xnyn  = H2dw(x, y, Q, Qx, Qy, Qxy, x0_xnyn, y0_xnyn);
	
	Q = (Q0_xpyp + Q0_xnyn + Q0_xpyn + Q0_xnyp)/4;
	Qx = (Q0_xpyp - Q0_xnyn + Q0_xpyn - Q0_xnyp)/(4*eh);
	Qy = (Q0_xpyp - Q0_xnyn - Q0_xpyn + Q0_xnyp)/(4*eh);
	Qxy = (Q0_xpyp + Q0_xnyn - Q0_xpyn - Q0_xnyp)/(4*eh*eh);
	
	t = t+dt;
	
	if(Ix > 129)
		t
	end;
	
	if(ifPlot)
		contourf(x, y, Q, 10);
		axis([xMin, xMax, yMin, yMax]);
		drawnow;
		print( strcat('Images/Image_',num2str(Ix), '_', num2str(ctr),'.png'));
		ctr = ctr+1;
	end;
		
end;
