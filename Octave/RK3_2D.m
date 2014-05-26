function [x, y] = RK3_2D(x0, y0, t0, dt, u, v)

x1 = x0 + dt*u(x0, y0, t0);
y1 = y0 + dt*v(x0, y0, t0);

x2 = x0 + dt*(0.25*u(x0, y0, t0) + 0.25*u(x1, y1, t0+dt));
y2 = y0 + dt*(0.25*v(x0, y0, t0) + 0.25*v(x1, y1, t0+dt));

x = x0 + dt*(u(x0, y0, t0)/6 + u(x1, y1, t0+dt)/6 + u(x2, y2, t0+0.5*dt)*2/3);
y = y0 + dt*(v(x0, y0, t0)/6 + v(x1, y1, t0+dt)/6 + v(x2, y2, t0+0.5*dt)*2/3);


