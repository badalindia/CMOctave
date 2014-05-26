function G = Diffeo2D(x, y, F, Ux, Uy, Nc, Nf, t0, tf, dt, ifPlot)

%evolves the function F 
%using velocity field Ux, Uy
%through diffeomorphism. 

page_output_immediately(true);

[Iy Ix] = size(x);
xMin = x(1, 1);
xMax = x(1, Ix);
xDelta = xMax - xMin;
yMin = y(1, 1);
yMax = y(Iy, 1);
yDelta = yMax - yMin;
h = x(1, 2) - x(1, 1);

%diffeomorphism grid

%coarse grid
Bc = linspace(xMin, xMax, Nc+1);
hc = Bc(2) - Bc(1);
[Xc, Yc] = [H, V] = meshgrid(Bc, Bc);
Hx = Vy = H.*0 + 1;
Hxy = Vxy = Hy = Vx = H.*0;

%fine grid
Bf = linspace(xMin, xMax, Nf+1);
hf = Bf(2) - Bf(1);
[Xf, Yf] = [Hf, Vf] = meshgrid(Bf, Bf);
Hxf = Vyf = Hf.*0 + 1;
Hxyf = Vxyf = Hyf = Vxf = Hf.*0;

t = t0;
i = j = [2:Nc];
G = F(x, y);
ctr = 0;

%some control parameters
nC = 10;
gridShiftN = 10;

if(ifPlot)
	contourf(x, y, G, nC);
	axis([xMin, xMax, yMin, yMax]);
	shading flat;
	drawnow;
	print( strcat('Images/Diffeo_',num2str(Nc), '_',num2str(Nf), '_', num2str(ctr),'.png'));
end;

ctr = ctr+1;

while t<tf;

	if(t+dt>tf)
		dt = tf - t;
	end;
	
	%backward integration
	eh = 10^-6;
	[x0_xpyp, y0_xpyp] = RK3_2D(Xc+eh, Yc+eh, t+dt, -dt, Ux, Uy);
	[x0_xpyn, y0_xpyn] = RK3_2D(Xc+eh, Yc-eh, t+dt, -dt, Ux, Uy);
	[x0_xnyp, y0_xnyp] = RK3_2D(Xc-eh, Yc+eh, t+dt, -dt, Ux, Uy);
	[x0_xnyn, y0_xnyn] = RK3_2D(Xc-eh, Yc-eh, t+dt, -dt, Ux, Uy);
	
	H0_xpyp  = H2dw(Xc, Yc, H, Hx, Hy, Hxy, x0_xpyp, y0_xpyp);
	H0_xpyn  = H2dw(Xc, Yc, H, Hx, Hy, Hxy, x0_xpyn, y0_xpyn);
	H0_xnyp  = H2dw(Xc, Yc, H, Hx, Hy, Hxy, x0_xnyp, y0_xnyp);
	H0_xnyn  = H2dw(Xc, Yc, H, Hx, Hy, Hxy, x0_xnyn, y0_xnyn);
	
	V0_xpyp  = H2dw(Xc, Yc, V, Vx, Vy, Vxy, x0_xpyp, y0_xpyp);
	V0_xpyn  = H2dw(Xc, Yc, V, Vx, Vy, Vxy, x0_xpyn, y0_xpyn);
	V0_xnyp  = H2dw(Xc, Yc, V, Vx, Vy, Vxy, x0_xnyp, y0_xnyp);
	V0_xnyn  = H2dw(Xc, Yc, V, Vx, Vy, Vxy, x0_xnyn, y0_xnyn);
	
	H = (H0_xpyp + H0_xnyn + H0_xpyn + H0_xnyp)/4;
	Hx = (H0_xpyp - H0_xnyn + H0_xpyn - H0_xnyp)/(4*eh);
	Hy = (H0_xpyp - H0_xnyn - H0_xpyn + H0_xnyp)/(4*eh);
	Hxy = (H0_xpyp + H0_xnyn - H0_xpyn - H0_xnyp)/(4*eh*eh);
	
	V = (V0_xpyp + V0_xnyn + V0_xpyn + V0_xnyp)/4;
	Vx = (V0_xpyp - V0_xnyn + V0_xpyn - V0_xnyp)/(4*eh);
	Vy = (V0_xpyp - V0_xnyn - V0_xpyn + V0_xnyp)/(4*eh);
	Vxy = (V0_xpyp + V0_xnyn - V0_xpyn - V0_xnyp)/(4*eh*eh);
	
	t = t+dt;
	
	#{
	%plotting
	if(ifPlot)
	
		%advected field
		_x = H2dw(Xc, Yc, H, Hx, Hy, Hxy, x, y);
		_y = H2dw(Xc, Yc, V, Vx, Vy, Vxy, x, y);
		
		__x = H2dw(Xf, Yf, Hf, Hxf, Hyf, Hxyf, _x, _y);
		__y = H2dw(Xf, Yf, Vf, Vxf, Vyf, Vxyf, _x, _y);
		
			G = F(__x, __y);

	
		contourf(x, y, G, nC);
		axis([xMin, xMax, yMin, yMax]);
		shading flat;
		drawnow;
		print( strcat('Images/Diffeo_',num2str(Nc), '_',num2str(Nf), '_', num2str(ctr),'.png'));
	end;
	#}
	
	if(mod(ctr,gridShiftN)==0)
	
	
		[Hf_t, Hxf_t, Hyf_t, Hxyf_t, Vf_t, Vxf_t, Vyf_t, Vxyf_t] = composeMaps(	Xf, Yf, xDelta, yDelta,
																Xc, Yc, H, Hx, Hy, Hxy, V, Vx, Vy, Vxy, 
																Xf, Yf, Hf, Hxf, Hyf, Hxyf, Vf, Vxf, Vyf, Vxyf);
	
		Hf = Hf_t;
		Hxf = Hxf_t;
		Hyf = Hyf_t;
		Hxyf = Hxyf_t;
		Vf = Vf_t;
		Vxf = Vxf_t;
		Vyf = Vyf_t;
		Vxyf = Vxyf_t;
		
		[Xc, Yc] = [H, V] = meshgrid(Bc, Bc);
		Hx = Vy = H.*0 + 1;
		Hxy = Vxy = Hy = Vx = H.*0;
		t
	end;
	
	ctr = ctr+1;
end


_x = H2dw(Xc, Yc, H, Hx, Hy, Hxy, x, y);
_y = H2dw(Xc, Yc, V, Vx, Vy, Vxy, x, y);
		
__x = H2dw(Xf, Yf, Hf, Hxf, Hyf, Hxyf, _x, _y);
__y = H2dw(Xf, Yf, Vf, Vxf, Vyf, Vxyf, _x, _y);

	G = F(__x, __y);
