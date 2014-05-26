function [xn, yn] = fpi(Xc, Yc, x, y, t, dt, A, Ax, Ay, Axy)

xn = x;
yn = y;

n = 1;
ctr = 0;
while e>10^-15 && ctr<20
[Axn, Ayn] = H2dw_dxdy(Xc, Yc, A, Ax, Ay, Axy, xn, yn);
Ux = -Ayn;
Uy = Axn;

e = max(max(max(abs(xn + dt.*Ux - x))), max(max(abs(yn + dt.*Uy - y))));
ctr = ctr+1;

xn = x - dt.*Ux;
yn = y - dt.*Uy;
end

if(e>10^-5 && ctr==20)
	printf("Fixed Point Iteration not converging\n");
end

