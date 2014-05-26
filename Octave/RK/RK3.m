function ans = RK3(x0, t0, dt, v)

x1 = x0 + dt*v(x0, t0);
x2 = x0 + dt*(0.25*v(x0, t0) + 0.25*v(x1, t0+dt));
ans = x0 + dt*(v(x0,t0)/6 + v(x1, t0+dt)/6 + v(x2, t0 + 0.5*dt)*2/3);
