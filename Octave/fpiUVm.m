function [xn, yn] = fpiUVm(Xc, Yc, x, y, t, dt, U, Ux, Uy, Uxy, V, Vx, Vy, Vxy)

xn = x;
yn = y;

n = 1;
ctr = 0;
while e>10^-15 && ctr<5
Un = H2dwm(Xc, Yc, U, Ux, Uy, Uxy, xn, yn);
Vn = H2dwm(Xc, Yc, V, Vx, Vy, Vxy, xn, yn);

e = max(max(max(abs(xn + dt.*Un - x))), max(max(abs(yn + dt.*Vn - y))));
ctr = ctr+1;

xn = x - dt.*Un;
yn = y - dt.*Vn;
end

if(e>10^-5 && ctr==20)
	printf("Fixed Point Iteration not converging\n");
end

