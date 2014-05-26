Nc = 32;
Nf = 64;
Np = 128;

Xc = linspace(0, 1, Nc+1);
hc = Xc(2) - Xc(1);

Xf = linspace(0, 1, Nf+1);
Xs = linspace(0, 1, Np+1);

%1D diffeo
Hc = Xc;
Hxc = Xc.*0 + 1;

Hf = Xf;
Hxf = Xf.*0 + 1;


%initial conditions
Q  = @(x) sin(2*pi*x);				%scalar field

%prescribed velocity
U = @(X, T) sin(2*pi*X)*cos(T);			%velocity field

t0 = 0;
tf = 2*pi;
dt = tf/200;
t = t0;

ep = 10^-6;

while t<tf

	if(t+dt>tf)
		dt = tf - t;
	end;
	
	
	x0p = RK3(Xc+ep, t+dt, -dt, U);
	x0n = RK3(Xc-ep, t+dt, -dt, U);
	
	Hp = H1dw(Xc, Hc, Hxc, x0p);
	Hn = H1dw(Xc, Hc, Hxc, x0n);
	
	Hc = (Hp+Hn)/2;
	Hxc = (Hp-Hn)/(2*ep);
	
	
	t = t+dt;
	
	_x = H1dw(Xc, Hc, Hxc, Xs);
	_Q = Q(_x);
	
	plot(Xs, _Q, '-or');
	drawnow;
	
end;

