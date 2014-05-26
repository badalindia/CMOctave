n = [4:6];
N = 2.^n;
I = length(N);

%error
e = zeros(I, 1);
Dx = zeros(I, 1);

t0 = 0;
tf = pi;
Np = 128;

x = linspace(0, 1, Np+1);
[X, Y] = meshgrid(x);
%G0 = RHO(X, Y, 1);



Nc_ref = N(I)*2;
Nf_ref = N(I)*4;
G0 = advection(Nc_ref, Nf_ref, Np, t0, tf, pi/(8*Nc_ref));
%G0 = dlmread("Images/rho_128_0.blk", "\n");

for i=[1:I]

	x = linspace(0, 1, N(i)+1);
	h = x(2) - x(1);

	G = advection(N(i), N(i)*2, Np, t0, tf, pi/(8*N(i)));

	e(i) = max(max(abs(G - G0)));
	Dx(i) = h;

end;

e
Dx
plot(log(Dx), log(e), '-o');
Accuracy = polyfit(log(Dx), log(e), 1)(1)
