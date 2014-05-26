%The following code shows how a 2D Hermite interpolant can be modified 
% to attain monotonicity

page_output_immediately(true);

x0 = 0; 
x1 = 1;
h = x1 - x0;

N = 16;
x = linspace(x0, x1, N+1);
x = x./h;
[X, Y] = meshgrid(x);

Hf = @(x) (1 - 3.*x.^2 + 2.*x.^3);
Hg = @(x) x.*(1 - x).^2;

i = [1:(N+1)*(N+1)];
F = [	Hf(X(i)) .* Hf(Y(i));				Hf(X(i)) .* Hf(1-Y(i)); 				Hf(X(i)) .* Hg(Y(i)); 			Hf(X(i)) .* -Hg(1-Y(i));
		Hf(1-X(i)) .* Hf(Y(i)); 			Hf(1-X(i)) .* Hf(1-Y(i));				Hf(1-X(i)) .* Hg(Y(i));			Hf(1-X(i)) .* -Hg(1-Y(i));
		Hg(X(i)) .* Hf(Y(i)); 			    Hg(X(i)) .* Hf(1-Y(i)); 				Hg(X(i)) .* Hg(Y(i)); 			Hg(X(i)) .* -Hg(1-Y(i));
		-Hg(1-X(i)) .* Hf(Y(i)); 			-Hg(1-X(i)) .* Hf(1-Y(i)); 			-Hg(1-X(i)) .* Hg(Y(i)); 			-Hg(1-X(i)) .* -Hg(1-Y(i))
		];

y  = zeros(N+1, N+1);
y2 = zeros(N+1, N+1);

% loop for generating numerous random interpolants to assess overshooting
shootCtr1 = 0;
shootCtr2 = 0;
plotCtr = 0;
for ctr = [1:50]

%random initialization of Interpolant node values
f00 = rand(); f10 = rand(); f01 = rand(); f11 = rand(); 
fx00 = rand(); fx10 = rand(); fx01 = rand(); fx11 = rand(); 
fy00 = rand(); fy10 = rand(); fy01 = rand(); fy11 = rand(); 
fxy00 = rand(); fxy10 = rand(); fxy01 = rand(); fxy11 = rand(); 


M = [	f00;		f01;		h*fy00;		h*fy01;
		f10;		f11;		h*fy10;		h*fy11;
		h*fx00;		h*fx01;		h*h*fxy00;	h*h*fxy01; 
		h*fx10;		h*fx11;		h*h*fxy10;	h*h*fxy11
		];
y(i) = M'*F;

minF = min([f00, f10, f01, f11]);
maxF = max([f00, f10, f01, f11]);

overShoot1 = sum(sum(y>maxF)) + sum(sum(y<minF));
	
	
%correction for monotonicity
dX0 = f10 - f00;
dX1 = f11 - f01;
dY0 = f01 - f00;
dY1 = f11 - f10;
D = dX1 - dX0;

%step 1, making x interpolant along the two edges monotone
fx00 = Bound(fx00, 3/h*dX0);
fx10 = Bound(fx10, 3/h*dX0);
fx01 = Bound(fx01, 3/h*dX1);
fx11 = Bound(fx11, 3/h*dX1);


if(dY1*dY0 < 0)
	
	fy00 = fy01 = fy10 = fy11 = 0;
	fxy00 = fxy01 = fxy10 = fxy11 = 0;
	
elseif (dY1*dY0 >= 0)

	fy00 = Bound(fy00, 3/h*dY0);
	fy01 = Bound(fy01, 3/h*dY0);
	fy10 = Bound(fy10, 3/h*dY1);
	fy11 = Bound(fy11, 3/h*dY1);

	%stronger conditions
	[fx00, fx01] = Bound2(fx00, fx01, 3/h*dX0, 3/h*dX1);
	[fx10, fx11] = Bound2(fx10, fx11, 3/h*dX0, 3/h*dX1);
	[fy00, fy10] = Bound2(fy00, fy10, 3/h*dY0, 3/h*dY1);
	[fy01, fy11] = Bound2(fy01, fy11, 3/h*dY0, 3/h*dY1);
	
	%second order conditions
	fxy00 = Bound(fxy00, 3/h*(fy10 - fy00));
	fxy00 = Bound(fxy00, 3/h*(fx01 - fx00));
	
	fxy10 = Bound(fxy10, 3/h*(fy10 - fy00));
	fxy10 = Bound(fxy10, 3/h*(fx11 - fx10));
	
	fxy01 = Bound(fxy01, 3/h*(fy11 - fy01));
	fxy01 = Bound(fxy01, 3/h*(fx01 - fx00));
	
	fxy11 = Bound(fxy11, 3/h*(fy11 - fy01));
	fxy11 = Bound(fxy11, 3/h*(fx11 - fx10));
	
end;


M = [	f00;		f01;		h*fy00;		h*fy01;
		f10;		f11;		h*fy10;		h*fy11;
		h*fx00;		h*fx01;		h*h*fxy00;	h*h*fxy01; 
		h*fx10;		h*fx11;		h*h*fxy10;	h*h*fxy11
		];
y2(i) = M'*F;

overShoot2 = sum(sum(y2>maxF)) + sum(sum(y2<minF));


if(overShoot1 > 0)
	shootCtr1 = shootCtr1 + 1;
end;

if(overShoot2 > 0)
	shootCtr2 = shootCtr2 + 1;
end;


	subplot(1, 2, 1)
	contourf(X, Y, y, 20);
	axis([0, 1, 0, 1], "square");
	caxis([minF, maxF]);
	subplot(1, 2, 2)
	contourf(X, Y, y2, 20);
	axis([0, 1, 0, 1], "square");
	caxis([minF, maxF]);
	print(strcat('ModifiedHermiteContour2D_', num2str(plotCtr),'.png'));
	plotCtr = plotCtr+1;

end;

shootCtr1
shootCtr2
