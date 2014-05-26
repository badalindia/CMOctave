clear;

%defining testing function
C1 = 1;
C2 = 2;
FX = @(x,y) sin(C1.*x).*sin(C2.*y);
FXx = @(x,y) C1.*cos(C1.*x).*sin(C2.*y);
FXy = @(x,y) C2.*sin(C1.*x).*cos(C2.*y);
FXxy = @(x,y) C1.*C2.*cos(C1.*x).*cos(C2.*y);

FY = @(x,y) cos(C1.*x).*cos(C2.*y);
FYx = @(x,y) -C1.*sin(C1.*x).*cos(C2.*y);
FYy = @(x,y) -C2.*cos(C1.*x).*sin(C2.*y);
FYxy = @(x,y) C1.*C2.*sin(C1.*x).*sin(C2.*y);


D1 = 2;
D2 = 1;
FH = @(x,y) sin(D1.*x).*sin(D2.*y);
FHx = @(x,y) D1.*cos(D1.*x).*sin(D2.*y);
FHy = @(x,y) D2.*sin(D1.*x).*cos(D2.*y);
FHxy = @(x,y) D1.*D2.*cos(D1.*x).*cos(D2.*y);

FC = @(x, y) FH(FX(x,y), FY(x, y));
FCx = @(x, y) C1*(D1*cos(C1*x).*sin(C2*y).*cos(D1*sin(C1*x).*sin(C2*y)).*sin(D2*cos(C1*x).*cos(C2*y))-D2*sin(C1*x).*cos(C2*y).*sin(D1*sin(C1*x).*sin(C2*y)).*cos(D2*cos(C1*x).*cos(C2*y)));
FCy = @(x, y) C2*(D1*sin(C1*x).*cos(C2*y).*cos(D1*sin(C1*x).*sin(C2*y)).*sin(D2*cos(C1*x).*cos(C2*y))-D2*cos(C1*x).*sin(C2*y).*sin(D1*sin(C1*x).*sin(C2*y)).*cos(D2*cos(C1*x).*cos(C2*y)));
FCxy = @(x, y) (
-C1*C2*D1^2.*sin(C1*x).*cos(C1*x).*sin(C2*y).*cos(C2*y).*sin(D1*sin(C1*x).*sin(C2*y)).*sin(D2*cos(C1*x).*cos(C2*y))
-C1*C2*D2^2.*sin(C1*x).*cos(C1*x).*sin(C2*y).*cos(C2*y).*sin(D1*sin(C1*x).*sin(C2*y)).*sin(D2*cos(C1*x).*cos(C2*y))
-C1*C2*D1*D2*(sin(C1*x).*cos(C2*y)).^2.*cos(D2*cos(C1*x).*cos(C2*y)).*cos(D1*sin(C1*x).*sin(C2*y))
-C1*C2*D1*D2*(cos(C1*x).*sin(C2*y)).^2.*cos(D2*cos(C1*x).*cos(C2*y)).*cos(D1*sin(C1*x).*sin(C2*y))
+C1*C2*D1*cos(C1*x).*cos(C2*y).*cos(D1*sin(C1*x).*sin(C2*y)).*sin(D2*cos(C1*x).*cos(C2*y))
+C1*C2*D2*sin(C1*x).*sin(C2*y).*sin(D1*sin(C1*x).*sin(C2*y)).*cos(D2*cos(C1*x).*cos(C2*y)));

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

%error 
e = zeros(I, 1);
d = zeros(I, 1);

			
for i=[1:I]

	x = linspace(0, 2*pi, N(i)+1);
	[X, Y] = meshgrid(x);

	fX = FX(X, Y);
	fXx = FXx(X, Y);
	fXy = FXy(X, Y);
	fXxy = FXxy(X, Y);

	fY = FY(X, Y);
	fYx = FYx(X, Y);
	fYy = FYy(X, Y);
	fYxy = FYxy(X, Y);

	fH = FH(X, Y);
	fHx = FHx(X, Y);
	fHy = FHy(X, Y);
	fHxy = FHxy(X, Y);
	
		%random points
		XR = rand(Nr, Nr)*(xr1 - xr0) + xr0;
		YR = rand(Nr, Nr)*(xr1 - xr0) + xr0;

	[fC, fCx, fCy, fCxy, gC, gCx, gCy, gCxy] = composeMapsA(XR, YR, X, Y, fX, fXx, fXy, fXxy, fY, fYx, fYy, fYxy, X, Y, fH, fHx, fHy, fHxy, fH, fHx, fHy, fHxy);
	
	e(i) = max(max(abs(gCxy -  FCxy(XR, YR) )));
	d(i) = x(2) - x(1);

end;

plot(log(d), log(e), '-ro');
Accuracy = polyfit(log(d), log(e), 1)(1)



