%testing GALS 1D using exact solution

page_output_immediately(true);

n = [4];
N = 2.^n;
I = length(N);

%error
e = zeros(I, 1);
Dx = zeros(I, 1);

for i=[1:I]

	Np = N(i)
	x = linspace(0, 1, N(i)+1);
	h = x(2) - x(1);

	Q  = @(x) sin(2*pi*x + pi*h);				%scalar field
	U  = @(X, T) 0.2*sin(2*pi*X)*cos(T);			%velocity field
	
	t0 = 0;
	tf = 4*pi;
	
	_Q = Diffeo1D(N(i), x, Q, U, t0, tf, h/2, true);
	
	e(i) = max(abs(_Q - Q(x)));
	Dx(i) = h;

end

if I>1
	plot(log(Dx), log(e), '-or');
	Accuracy = polyfit(log(Dx), log(e), 1)(1)
end;
