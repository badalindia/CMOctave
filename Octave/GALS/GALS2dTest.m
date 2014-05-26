%testing GALS 2D using exact solution

n = [8];
N = 2.^n;
I = length(N);

%error
e = zeros(I, 1);
Dx = zeros(I, 1);

for i=[1:I]

	x = linspace(0, 1, N(i)+1);

	[X, Y] = meshgrid(x, x);

	C = 2;
	Q  = sin(C*pi*X).*sin(C*pi*Y);				
	Qx = C*pi*cos(C*pi*X).*sin(C*pi*Y);
	Qy = C*pi*sin(C*pi*X).*cos(C*pi*Y);
	Qxy= C*C*pi*pi*cos(C*pi*X).*cos(C*pi*Y);

	U = @(X, Y, T) -cos(T).*(sin(pi.*X).^2).*sin(2.*pi.*Y);
	V = @(X, Y, T)  cos(T).*(sin(pi.*Y).^2).*sin(2.*pi.*X);

	t0 = 0;
	tf = pi;
	h = x(2) - x(1);
	
	QExact = Q;

	Q = GALS2d(X, Y, Q, Qx, Qy, Qxy, U, V, t0, tf, h/2, true);
	
	e(i) = max(max(abs(Q - QExact)));
	Dx(i) = h;

end;

plot(log(Dx), log(e), '-o');
Accuracy = polyfit(log(Dx), log(e), 1)(1)
log(Dx)
log(e)
