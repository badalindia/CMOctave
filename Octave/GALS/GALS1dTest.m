%testing GALS 1D using exact solution

n = [6:10];
N = 2.^n;
I = length(N);

%error
e = zeros(I, 1);
Dx = zeros(I, 1);

for i=[1:I]

	x = linspace(-pi, pi, N(i)+1);
	
	Q  = sin(x);				%scalar field
	Qx = cos(x); 				%gradient
	U = @(X, T) sin(X)*cos(T);			%velocity field
	
	t0 = 0;
	tf = 2*pi;
	h = x(2) - x(1);

	QExact = sin(x);

	[Q Qx] = GALS1d(x, Q, Qx, U, t0, tf, h/2, false);
	
	e(i) = max(abs(Q - QExact));
	Dx(i) = h;

end

plot(log(Dx), log(e));
Accuracy = polyfit(log(Dx), log(e), 1)(1)
