function [f, fx, fy, fxy] = ModifyHermiteGrid(x, y, f, fx, fy, fxy)

%Modified Hermite node values at grid to make it monotone


[Ix, Iy] = size(x);
h = x(1, 2) - x(1, 1);

for j = [1:Iy-1]
for i = [1:Ix-1]

	dX0 = f(j, i+1) - f(j, i);
	dX1 = f(j+1, i+1) - f(j+1, i);
	dY0 = f(j+1, i) - f(j, i);
	dY1 = f(j+1, i+1) - f(j, i+1);
	D = dX1 - dX0;

	%step 1, making x interpolant along the two edges monotone
	fx(j, i) 		= Bound(fx(j, i), 3/h*dX0);
	fx(j, i+1) 		= Bound(fx(j, i+1), 3/h*dX0);
	fx(j+1, i) 		= Bound(fx(j+1, i), 3/h*dX1);
	fx(j+1, i+1) 	= Bound(fx(j+1, i+1), 3/h*dX1);
	
	if(dY1*dY0 < 0)
	
	fy(j, i) = fy(j+1, i) = fy(j, i+1) = fy(j+1, i+1) = 0;
	fxy(j, i) = fxy(j+1, i) = fxy(j, i+1) = fxy(j+1, i+1) = 0;
	
	elseif (dY1*dY0 >= 0)

		fy(j, i) = Bound(fy(j, i), 3/h*dY0);
		fy(j+1, i) = Bound(fy(j+1, i), 3/h*dY0);
		fy(j, i+1) = Bound(fy(j, i+1), 3/h*dY1);
		fy(j+1, i+1) = Bound(fy(j+1, i+1), 3/h*dY1);

		%stronger conditions
		[fx(j, i), fx(j+1, i)] = Bound2(fx(j, i), fx(j+1, i), 3/h*dX0, 3/h*dX1);
		[fx(j, i+1), fx(j+1, i+1)] = Bound2(fx(j, i+1), fx(j+1, i+1), 3/h*dX0, 3/h*dX1);
		[fy(j, i), fy(j, i+1)] = Bound2(fy(j, i), fy(j, i+1), 3/h*dY0, 3/h*dY1);
		[fy(j+1, i), fy(j+1, i+1)] = Bound2(fy(j+1, i), fy(j+1, i+1), 3/h*dY0, 3/h*dY1);
		
		%second order conditions
		fxy(j, i) = Bound(fxy(j, i), 3/h*(fy(j, i+1) - fy(j, i)));
		fxy(j, i) = Bound(fxy(j, i), 3/h*(fx(j+1, i) - fx(j, i)));
		
		fxy(j, i+1) = Bound(fxy(j, i+1), 3/h*(fy(j, i+1) - fy(j, i)));
		fxy(j, i+1) = Bound(fxy(j, i+1), 3/h*(fx(j+1, i+1) - fx(j, i+1)));
		
		fxy(j+1, i) = Bound(fxy(j+1, i), 3/h*(fy(j+1, i+1) - fy(j+1, i)));
		fxy(j+1, i) = Bound(fxy(j+1, i), 3/h*(fx(j+1, i) - fx(j, i)));
		
		fxy(j+1, i+1) = Bound(fxy(j+1, i+1), 3/h*(fy(j+1, i+1) - fy(j+1, i)));
		fxy(j+1, i+1) = Bound(fxy(j+1, i+1), 3/h*(fx(j+1, i+1) - fx(j, i+1)));
		
	end;
	

end;
end;

fx(:, Ix) = fx(:, 1);
fx(Iy, :) = fx(1, :);

fy(:, Ix) = fy(:, 1);
fy(Iy, :) = fy(1, :);

fxy(:, Ix) = fxy(:, 1);
fxy(Iy, :) = fxy(1, :);



