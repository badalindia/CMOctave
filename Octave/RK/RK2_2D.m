function [x, y] = RK2_2D(x0, y0, t0, dt, u, v)

kx1 = u(x0, y0, t0);
ky1 = v(x0, y0, t0);

kx2 = u(x0 + 0.5*dt*kx1, y0 + 0.5*dt*ky1, t0 + 0.5*dt);
ky2 = v(x0 + 0.5*dt*kx1, y0 + 0.5*dt*ky1, t0 + 0.5*dt);

x = x0 + dt*kx2;
y = y0 + dt*ky2;
