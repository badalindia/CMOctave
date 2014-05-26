%The following code shows how a 1D Hermite interpolant can be modified 
% to attain monotonicity


x0 = 0;
x1 = rand();
h = x1 - x0;

N = 16;
x = linspace(x0, x1, N+1);
x = x./h;

f0 = rand();
f1 = rand();
fx0 = rand();
fx1 = rand();

Hf = @(x) (1 - 3.*x.^2 + 2.*x.^3);
Hg = @(x) x.*(1 - x).^2;

y = f0.*Hf(x) + f1.*Hf(1 - x) + h*(fx0.*Hg(x) - fx1.*Hg(1 - x));
overshoot1 = sum(y>max([f0, f1])) + sum(y<min([f0, f1]))

%Modification
delta = (f1 - f0);

if( delta >= 0)

	if(fx0 < 0)
		fx0 = 0;
	end;
	if(fx0 > 3*delta/h)
		fx0 = 2.9*delta/h;
	end;
	
	if(fx1 < 0)
		fx1 = 0;
	end;
	if(fx1 > 3*delta/h)
		fx1 = 2.9*delta/h;
	end;
	
end;

if( delta < 0)

	if(fx0 > 0)
		fx0 = 0;
	end;
	if(fx0 < 3*delta/h)
		fx0 = 2.99*delta/h;
	end;
	
	if(fx1 > 0)
		fx1 = 0;
	end;
	if(fx1 < 3*delta/h)
		fx1 = 2.99*delta/h;
	end;

end;

y2 = f0.*Hf(x) + f1.*Hf(1 - x) + h*(fx0.*Hg(x) - fx1.*Hg(1 - x));
overshoot2 = sum(y2>max([f0, f1])) + sum(y2<min([f0, f1]))


plot(x, y, '-or', x, y2, '-ob');






