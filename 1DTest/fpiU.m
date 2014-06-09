function xn = fpiU(Xc, x, t, dt, U, Ux)

xn = x;

n = 1;
ctr = 0;
while e>10^-15 && ctr<5
Un = H1dw(Xc, U, Ux, xn);

e = max(abs(xn + dt.*Un - x));
ctr = ctr+1;

xn = x - dt.*Un;
end

if(e>10^-5 && ctr==20)
	printf("Fixed Point Iteration not converging\n");
end
