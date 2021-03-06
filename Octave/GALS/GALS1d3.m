%1-d advection (without exact solution)

K = [4];
E = K.*0;
H = K.*0;
ctr = 1;

for k=K

	%number of grid points
	N = 2^k + 1;

	%initializing scalar field
	x = linspace(-pi, pi, N);
	h = x(2) - x(1);

	t = 0;
	tEnd = 100;

	Q  = sin(x);				%scalar field
	Qx = cos(x); 	%gradient of scalar field (remains constant)
	U = @(X, T) sin(X).*cos(T);			%velocity field

	%time stepping
	while t<tEnd;

		u = U(x, t);
		uMax = max(u);
		dt = h/(2*uMax);
		
		if(t+dt>100)
			dt = 100 - t;
		end;

		%backward integration
		x0 = RK3(x, t+dt, -dt, U);

		%trick to find derivative
		ex = 0.5*(sqrt(h))^3;
		x0p = RK3(x+ex, t+dt, -dt, U);
		x0n = RK3(x-ex, t+dt, -dt, U);

		[Q0, Qx0] = H1d(x, Q, Qx, x0);
		
		Q0p = H1d(x, Q, Qx, x0p);
		Q0n = H1d(x, Q, Qx, x0n);

		Q = Q0;
		Qx = (Q0p - Q0n)/(2*ex);

		t = t+dt;

		if N==65	%plot the first simulation
			plot(x, Q, '-r');
			axis([-pi, pi, 0, 2]);
			drawnow;
		end;	

	end;
	
	if k==K(1)
		Qr = Q;	%reference Q at highest grid resolution
	end;

	i1 = [1:N];
	i0 = [1:2^(K(1) - k):2^K(1) + 1];
	E(ctr) = max(abs(Q(i1) - Qr(i0)));
	H(ctr) = h;
	ctr = ctr+1;

end;

plot(log(H), log(E), '-o', log(H), 2*log(H), '-r', log(H), 3*log(H), '-r', log(H), 4*log(H), '-r');


