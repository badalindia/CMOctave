function G = advection(Nc, Nf, Np, t0, tf, dt)


Nc
Nf
Np

page_output_immediately(true);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Grid
dL = 1;
x = linspace(0, dL, Nc+1);
h = x(2) - x(1);
[X, Y] = meshgrid(x);

%sampling grid
[Xs, Ys] = meshgrid(linspace(0, dL, Np+1));


        %diffeomorphism grids

        %coarse grid
        Bc = linspace(0, dL, Nc+1);
        hc = Bc(2) - Bc(1);
        [Xc, Yc] = [Hc, Vc] = meshgrid(Bc, Bc);
        Hxc = Vyc = Hc.*0 + 1;
        Hxyc = Vxyc = Hyc = Vxc = Hc.*0;

        %fine grid
        Bf = linspace(0, dL, Nf+1);
        hf = Bf(2) - Bf(1);
        [Xf, Yf] = [Hf, Vf] = meshgrid(Bf, Bf);
        Hxf = Vyf = Hf.*0 + 1;
        Hxyf = Vxyf = Hyf = Vxf = Hf.*0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

coeff = 1.5;
eh = 10^-4;
gridShiftN = 10;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plotting parameters
plotCtr = 64;
ifPlot = true;
target = "Advection";
Res = '-r150';
nC = 50;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


t = t0;
ctr = 1;
ctr2 = 0;

UU = @(X, Y, T) -cos(T).*(sin(pi.*X).^2).*sin(2.*pi.*Y);
VV = @(X, Y, T)  cos(T).*(sin(pi.*Y).^2).*sin(2.*pi.*X);


while t < tf - 0.1*dt;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[x0_xpyp, y0_xpyp] = RK3_2D(Xc+eh, Yc+eh, t+dt, -dt, UU, VV);
[x0_xpyn, y0_xpyn] = RK3_2D(Xc+eh, Yc-eh, t+dt, -dt, UU, VV);
[x0_xnyp, y0_xnyp] = RK3_2D(Xc-eh, Yc+eh, t+dt, -dt, UU, VV);
[x0_xnyn, y0_xnyn] = RK3_2D(Xc-eh, Yc-eh, t+dt, -dt, UU, VV);

H0_xpyp  = H2dw(Xc, Yc, Hc, Hxc, Hyc, Hxyc, x0_xpyp, y0_xpyp);
H0_xpyn  = H2dw(Xc, Yc, Hc, Hxc, Hyc, Hxyc, x0_xpyn, y0_xpyn);
H0_xnyp  = H2dw(Xc, Yc, Hc, Hxc, Hyc, Hxyc, x0_xnyp, y0_xnyp);
H0_xnyn  = H2dw(Xc, Yc, Hc, Hxc, Hyc, Hxyc, x0_xnyn, y0_xnyn);

V0_xpyp  = H2dw(Xc, Yc, Vc, Vxc, Vyc, Vxyc, x0_xpyp, y0_xpyp);
V0_xpyn  = H2dw(Xc, Yc, Vc, Vxc, Vyc, Vxyc, x0_xpyn, y0_xpyn);
V0_xnyp  = H2dw(Xc, Yc, Vc, Vxc, Vyc, Vxyc, x0_xnyp, y0_xnyp);
V0_xnyn  = H2dw(Xc, Yc, Vc, Vxc, Vyc, Vxyc, x0_xnyn, y0_xnyn);

	
Hc = (H0_xpyp + H0_xnyn + H0_xpyn + H0_xnyp)/4;
Hxc = (H0_xpyp - H0_xnyn + H0_xpyn - H0_xnyp)/(4*eh);
Hyc = (H0_xpyp - H0_xnyn - H0_xpyn + H0_xnyp)/(4*eh);
Hxyc = (H0_xpyp + H0_xnyn - H0_xpyn - H0_xnyp)/(4*eh*eh);

Vc = (V0_xpyp + V0_xnyn + V0_xpyn + V0_xnyp)/4;
Vxc = (V0_xpyp - V0_xnyn + V0_xpyn - V0_xnyp)/(4*eh);
Vyc = (V0_xpyp - V0_xnyn - V0_xpyn + V0_xnyp)/(4*eh);
Vxyc = (V0_xpyp + V0_xnyn - V0_xpyn - V0_xnyp)/(4*eh*eh);

t = t+dt;
	
	%grid composition and resetting
	if(mod(ctr,gridShiftN)==0)
	[Hf, Hxf, Hyf, Hxyf, Vf, Vxf, Vyf, Vxyf] = composeMaps(Xf, Yf, 0, 0,
												Xc, Yc, Hc, Hxc, Hyc, Hxyc, Vc, Vxc, Vyc, Vxyc, 
												Xf, Yf, Hf, Hxf, Hyf, Hxyf, Vf, Vxf, Vyf, Vxyf);
												
	
	[Xc, Yc] = [Hc, Vc] = meshgrid(Bc, Bc);
	Hxc = Vyc = Hc.*0 + 1;
	Hxyc = Vxyc = Hyc = Vxc = Hc.*0;                
	t
	end;
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



if(mod(ctr, plotCtr)==0 && ifPlot)
        
        t

		%using sampling grid
		_x = H2dw(Xc, Yc, Hc, Hxc, Hyc, Hxyc, Xs, Ys);
		_y = H2dw(Xc, Yc, Vc, Vxc, Vyc, Vxyc, Xs, Ys);
		
		__x = H2dw(Xf, Yf, Hf, Hxf, Hyf, Hxyf, _x, _y);
		__y = H2dw(Xf, Yf, Vf, Vxf, Vyf, Vxyf, _x, _y);
		
		G = RHO(__x, __y, coeff);
			
			maxG = max(max(G));
			minG = min(min(G));
			
		
        hold off;
		contour(Xs, Ys, 2-G, [1,1], "linewidth", 2);
		axis([0, 1, 0, 1],"square");
		
		hold on;
		plot(Xs, Xs', 'r', "linewidth", 1);
		plot(Ys, Ys', 'r', "linewidth", 1);
		hold off;
		
		print( strcat('Images/rho/RHO_', num2str(Nc), '_', num2str(ctr2),'.png') , '-r10000');
		
		
		%dlmwrite(strcat('Images/rho_', num2str(Nc), '_', num2str(ctr2),'.blk'), 2 - G, "delimiter", "\n");
		
		#{
		%plotting transformed grid
        x1 = linspace(0, 1, Nc);
        x2 = linspace(0, 1, 256);
        
        	[Xv, Yv] = meshgrid(x1, x2);
			[Xh, Yh] = meshgrid(x2, x1);
			
			_Xv = H2dw(Xc, Yc, Hc, Hxc, Hyc, Hxyc, Xv, Yv);
			_Yv = H2dw(Xc, Yc, Vc, Vxc, Vyc, Vxyc, Xv, Yv);
			
			__Xv = H2dw(Xf, Yf, Hf, Hxf, Hyf, Hxyf, _Xv, _Yv);
			__Yv = H2dw(Xf, Yf, Vf, Vxf, Vyf, Vxyf, _Xv, _Yv);
			
			_Xh = H2dw(Xc, Yc, Hc, Hxc, Hyc, Hxyc, Xh, Yh);
			_Yh = H2dw(Xc, Yc, Vc, Vxc, Vyc, Vxyc, Xh, Yh);
			
			__Xh = H2dw(Xf, Yf, Hf, Hxf, Hyf, Hxyf, _Xh, _Yh);
			__Yh = H2dw(Xf, Yf, Vf, Vxf, Vyf, Vxyf, _Xh, _Yh);
			
			plot(__Xv, __Yv, 'r', "linewidth", 1);
			axis([0, 1, 0, 1], 'square');
			hold on;
			plot(__Xh', __Yh', 'r', "linewidth", 1);
			print( strcat('Images/Diffeo_', num2str(Nc), '_', num2str(ctr2),'.png') , '-r100');
			hold off;
		#}
				
			
		ctr2 = ctr2+1;	
        
end

ctr = ctr+1;

end;

ctr
t

#{
%using sampling grid
_x = H2dw(Xc, Yc, Hc, Hxc, Hyc, Hxyc, Xs, Ys);
_y = H2dw(Xc, Yc, Vc, Vxc, Vyc, Vxyc, Xs, Ys);

__x = H2dw(Xf, Yf, Hf, Hxf, Hyf, Hxyf, _x, _y);
__y = H2dw(Xf, Yf, Vf, Vxf, Vyf, Vxyf, _x, _y);
		
G = RHO(__x, __y, coeff);

dlmwrite(strcat('Images/rho_', num2str(Nc), '_', num2str(Nf), '_', num2str(Np), '_', num2str(ctr2),'.blk'), 2 - G, "delimiter", "\n");

csvwrite(strcat('Images/rho_', num2str(Nc), '_', num2str(Nf), '_', num2str(Np), '_', num2str(ctr2),'.csv'), 2 - G);
#}



