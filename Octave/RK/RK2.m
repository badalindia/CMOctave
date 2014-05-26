function ans = RK2(x0, t0, dt, v)

k1 = v(x0, t0);
k2 = v(x0 + 0.5*dt*k1, t0 + 0.5*dt);
ans = x0 + dt*k2;
