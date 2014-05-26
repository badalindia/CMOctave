%tests the order of accuracy of Hermite Interpolation formula


%testing function
C = 2;
f = @(x) sin(C*x);
fx = @(x) C*cos(C*x);

%test conditions
x0 = -pi;	
x1 = pi;
n = [3:12];		%time domain will divide into 2^n equal step sizes

N = 2.^n;
Dx = (x1 - x0)./N;
I = length(n);

ep = 0.001;
_x = linspace(x0+ep, x1-ep, 2^(n(I)+4));
fExact = f(_x);

%error 
e = zeros(I, 1);

for i=[1:I]

	%points of sampling
	xS = linspace(x0, x1, N(i)+1);	
	fS = f(xS);
	fxS = fx(xS);
	
	%hermite interpolation	
	h = H1d(xS, fS, fxS, _x);
	
	e(i) = max(abs(h - fExact));
	
end

plot(log(Dx), log(e), '-o');
Accuracy = polyfit(log(Dx'), log(e), 1)(1)
