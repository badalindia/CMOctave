
n = [3:6];
I = length(n);

%error
e = zeros(I, 1);
Dx = zeros(I, 1);

U = @(X, Y, T) -cos(T).*(sin(pi.*X).^2).*sin(2.*pi.*Y);
V = @(X, Y, T)  cos(T).*(sin(pi.*Y).^2).*sin(2.*pi.*X);

for i = [1:I]

	i
	Nc = 2^n(i)
	Nf = 2^(n(i)+1)
	Np = 2^(n(i)+1);

	x = linspace(0, 1, Np+1);

	[X, Y] = meshgrid(x, x);

	C = 2;
	Q = @(X, Y) sin(C*pi*X).*sin(C*pi*Y);	

	h = x(2) - x(1);

	Qh = Diffeo2D(X, Y, Q, U, V, Nc, Nf, 0, pi, h/2, false);
	Q0 = Q(X, Y);
	
	e(i) = max(max(abs(Qh - Q0)));
	Dx(i) = h;

end;

plot(log(Dx), log(e), '-o');
Accuracy = polyfit(log(Dx), log(e), 1)(1)
e
