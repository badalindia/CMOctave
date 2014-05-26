function ans = RK4(x0, t0, dt, v)

k1 = v(x0, t0);
k2 = v(x0 + 0.5*dt*k1, t0 + 0.5*dt);
k3 = v(x0 + 0.5*dt*k2, t0 + 0.5*dt);
k4 = v(x0 + dt*k3, t0 + dt);
ans = x0 + dt*(k1 + 2*k2 + 2*k3 + k4)/6;

