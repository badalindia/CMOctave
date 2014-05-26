%tests the order of accuracy of 2-D Hermite Interpolation formula
clear;

%testing function
C1 = 2;
C2 = 3;
C3 = 1.5;
C4 = 2.1;
f = @(x,y) sin(C1.*x).*sin(C2.*y) + C1*x + C2*y;
fx = @(x,y) C1.*cos(C1.*x).*sin(C2.*y) + C1;
fy = @(x,y) C2.*sin(C1.*x).*cos(C2.*y) + C2;
fxy = @(x,y) C1.*C2.*cos(C1.*x).*cos(C2.*y);

%f = @(x,y) x;
%fx = @(x,y) x.*0 + 1;
%fy = @(x,y) x.*0;
%fxy = @(x,y) x.*0;

%test conditions
n = [4:9];		%time domain will divide into 2^n equal step sizes
N = 2.^n;
I = length(n);
Nr = 10;	%count of random numbers

%grid range
x0 = 0;
x1 = 2*pi;

%random points range
xr0 = -2*pi;
xr1 = 4*pi;

_xR = linspace(xr0, xr1, 3*N(1)+1); 
[_XR, _YR] = meshgrid(_xR);


%error 
e = zeros(I, 1);
d = zeros(I, 1);

for i=[1:I]

	%periodic grid
	xS = linspace(x0, x1, N(i)+1); 
	h = xS(2) - xS(1);
	[XS, YS] = meshgrid(xS, xS);
	
	%function values at periodic grid
	fS = f(XS, YS);
	fxS = fx(XS, YS);
	fyS = fy(XS, YS);
	fxyS = fxy(XS, YS);
	
	%random points
	XR = rand(Nr, Nr)*(xr1 - xr0) + xr0;
	YR = rand(Nr, Nr)*(xr1 - xr0) + xr0;
	
	
	#{
	plot(_XR, _XR', '-g');
	axis([min(x0, xr0), max(x1, xr1), min(x0, xr0), max(x1, xr1)]);
	hold on;
	plot(_YR, _YR', '-g');
	plot(XS, XS', '-r');
	plot(YS, YS', '-r');
	scatter(XR, YR);
	hold off;
	drawnow;
	#}
	
	
	%hermite interpolation	
	[fh fhx fhy fhxy] = H2dwF(XS, YS, fS, fxS, fyS, fxyS, XR, YR);
	fxyExact = fxy(XR, YR);	
	
	e(i) = max(max(abs(fhxy - fxyExact)));
	d(i) = h;
	
	
end

plot(log(d), log(e), '-ro');
Accuracy = polyfit(log(d), log(e), 1)(1)

