%testing GALS 1D using exact solution

page_output_immediately(true);

n = [5:8];
N = 2.^n;
I = length(N);

%error
e = zeros(I, 1);
Dx = zeros(I, 1);

for i=[1:I]

	Np = N(i)
	x = linspace(0, 1, N(i)+1);
	
	Q  = @(x) sin(2*pi*x);				%scalar field
	U  = @(X, T) 0.2*sin(2*pi*X)*cos(T);			%velocity field
	
	t0 = 0;
	tf = 2*pi;
	h = x(2) - x(1);

	_Q = Diffeo1D(N(i), x, Q, U, t0, tf, h/2, false);
	
	e(i) = max(abs(_Q - Q(x)));
	Dx(i) = h;

end

if I>1
	plot(log(Dx), log(e), '-or');
	Accuracy = polyfit(log(Dx), log(e), 1)(1)
end;

I0yCell= mod(I0yPeriod,2)*(-(mod(I0y-1,L-1)+1)-(I0y-sign(I0yPeriod)*(L-1)*(I0yPeriod-1))+((sign(I0yPeriod)+1)*7+2)-)  +  mod(I0y-1,L-1)+1;
