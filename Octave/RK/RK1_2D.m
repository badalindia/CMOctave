function [x, y] = RK1_2D(x0, y0, t0, dt, u, v)

kx1 = u(x0, y0, t0);
ky1 = v(x0, y0, t0);

x = x0 + dt*kx1;
y = y0 + dt*ky1;
