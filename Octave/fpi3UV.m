function [xf2, yf2] = fpi3UV(Xc, Yc, x, y, t, dt, U, Ux, Uy, Uxy, V, Vx, Vy, Vxy, U_p, Ux_p, Uy_p, Uxy_p, V_p, Vx_p, Vy_p, Vxy_p)

%foot point at t-dt
xf1 = x;
yf1 = y;

%foot point at t
xf2 = x;
yf2 = y;

e = 1;
ctr = 0;
while e>10^-15 && ctr<20

%velocity at points at t-dt
U1 = H2dw(Xc, Yc, U_p, Ux_p, Uy_p, Uxy_p, xf1, yf1);
V1 = H2dw(Xc, Yc, V_p, Vx_p, Vy_p, Vxy_p, xf1, yf1);

xf1 = xf2 - dt*U1;
yf1 = yf2 - dt*V1;

%velocity at points at t-dt
U1 = H2dw(Xc, Yc, U_p, Ux_p, Uy_p, Uxy_p, xf1, yf1);
V1 = H2dw(Xc, Yc, V_p, Vx_p, Vy_p, Vxy_p, xf1, yf1);

%velocity at points at t-dt
U2 = H2dw(Xc, Yc, U, Ux, Uy, Uxy, xf2, yf2);
V2 = H2dw(Xc, Yc, V, Vx, Vy, Vxy, xf2, yf2);

xf2 = x - dt*(1.5*U2 - 0.5*U1);
yf2 = y - dt*(1.5*V2 - 0.5*V1);


e = max(max(max(abs(x - dt*(1.5*U2 - 0.5*U1) - xf2))), max(max(abs(y - dt*(1.5*V2 - 0.5*V1) - yf2))));
ctr = ctr+1;

end;

if(e>10^-5 && ctr==20)
	printf("Fixed Point Iteration not converging\n");
end
