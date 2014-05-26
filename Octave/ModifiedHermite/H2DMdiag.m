%The following code shows how a 2D Hermite interpolant can be modified 
% to attain monotonicity

%tests only along diagonals

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



% loop for generating numerous random interpolants to assess overshooting
shootCtr1 = 0;
shootCtr2 = 0;
plotCtr = 0;
for ctr = [1:500]


%random initialization of Interpolant node values
f00 = rand(); f10 = rand(); f01 = rand(); f11 = rand(); 
fx00 = rand(); fx10 = rand(); fx01 = rand(); fx11 = rand(); 
fy00 = rand(); fy10 = rand(); fy01 = rand(); fy11 = rand(); 
fxy00 = rand(); fxy10 = rand(); fxy01 = rand(); fxy11 = rand(); 


y = ( f00*Hf(x).* Hf(x) + f10*Hf(1-x).* Hf(x) + f01*Hf(x).* Hf(1-x) + f11*Hf(1-x).* Hf(1-x)
	+ h*(fx00*Hg(x).* Hf(x) + fx10*Hg(1-x).* Hf(x) + fx01*Hg(x).* Hf(1-x) + fx11*Hg(1-x).* Hf(1-x))
	+ h*(fy00*Hf(x).* Hg(x) + fy10*Hf(1-x).* Hg(x) + fy01*Hf(x).* Hg(1-x) + fy11*Hf(1-x).* Hg(1-x))
	+ h*h*(fxy00*Hg(x).* Hg(x) + fxy10*Hg(1-x).* Hg(x) + fxy01*Hg(x).* Hg(1-x) + fxy11*Hg(1-x).* Hg(1-x))
	);
overShoot1 = sum(y>max([f00, f10, f01, f11])) + sum(y<min([f00, f10, f01, f11]));
	
	
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


y2 = ( f00*Hf(x).* Hf(x) + f10*Hf(1-x).* Hf(x) + f01*Hf(x).* Hf(1-x) + f11*Hf(1-x).* Hf(1-x)
	+ h*(fx00*Hg(x).* Hf(x) - fx10*Hg(1-x).* Hf(x) + fx01*Hg(x).* Hf(1-x) - fx11*Hg(1-x).* Hf(1-x))
	+ h*(fy00*Hf(x).* Hg(x) + fy10*Hf(1-x).* Hg(x) - fy01*Hf(x).* Hg(1-x) - fy11*Hf(1-x).* Hg(1-x))
	+ h*h*(fxy00*Hg(x).* Hg(x) - fxy10*Hg(1-x).* Hg(x) - fxy01*Hg(x).* Hg(1-x) + fxy11*Hg(1-x).* Hg(1-x))
	);

overShoot2 = sum(y2>max([f00, f10, f01, f11])) + sum(y2<min([f00, f10, f01, f11]));


if(overShoot1 > 0)
	shootCtr1 = shootCtr1 + 1;
end;

if(overShoot2 > 0)
	shootCtr2 = shootCtr2 + 1;
end;


if(mod(ctr, 500) == 0)
	ctr
	shootCtr1
	shootCtr2
end;

	

#{
%print the image if the overshooting is dramatic
if(overShoot1 > 5)
	plot(x, y, '-or', x, y2, '-ob', x, x.*0 + max([f00, f10, f01, f11]), '-g', x, x.*0 + min([f00, f10, f01, f11]), '-g');
	print(strcat('ModifiedHermite2D_', num2str(plotCtr),'.png'));
	plotCtr = plotCtr+1;
end;
#}

end;

shootCtr1
shootCtr2
