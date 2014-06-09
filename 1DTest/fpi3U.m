function xf = fpi3U(Xc, x, t, dt, U, Ux, U_p, Ux_p)

%foot point at t-dt
xf1 = x;

%foot point at t
xf2 = x;

e = 1;
ctr = 0;
while e>10^-15 && ctr<20

%velocity at points at t-dt
U1 = H1dw(Xc, U_p, Ux_p, xf1);

xf1 = xf2 - dt*U1;

%velocity at points at t-dt
U1 = H1dw(Xc, U_p, Ux_p, xf1);

%velocity at points at t-dt
U2 = H1dw(Xc, U, Ux, xf2);

xf2 = x - dt*(1.5*U2 - 0.5*U1);

e = max(abs(x - dt*(1.5*U2 - 0.5*U1) - xf2));
ctr = ctr+1;

end;

if(e>10^-5 && ctr==20)
	printf("Fixed Point Iteration not converging\n");
end

xf = xf2;
