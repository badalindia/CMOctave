clear;

page_output_immediately(true);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Initial conditions
Nc = 8;
Nf = 16;
Np = 32;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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

t = 0;
dt = 0.01;
tf = 501*dt;
ctr = 1;
eh = 10^-4;
gridShiftN = 10;


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

	[Hc, Hxc, Hyc, Hxyc] = ModifyHermiteGrid(Xc, Yc, Hc, Hxc, Hyc, Hxyc);
	[Vc, Vxc, Vyc, Vxyc] = ModifyHermiteGrid(Xc, Yc, Vc, Vxc, Vyc, Vxyc);
	

t = t+dt;
	
	%grid composition and resetting
	if(mod(ctr,gridShiftN)==0)
	[Hf, Hxf, Hyf, Hxyf, Vf, Vxf, Vyf, Vxyf] = composeMaps(Xf, Yf, dL, dL,
												Xc, Yc, Hc, Hxc, Hyc, Hxyc, Vc, Vxc, Vyc, Vxyc, 
												Xf, Yf, Hf, Hxf, Hyf, Hxyf, Vf, Vxf, Vyf, Vxyf);

	[Hf, Hxf, Hyf, Hxyf] = ModifyHermiteGrid(Xf, Yf, Hf, Hxf, Hyf, Hxyf);
	[Vf, Vxf, Vyf, Vxyf] = ModifyHermiteGrid(Xf, Yf, Vf, Vxf, Vyf, Vxyf);
	
	[Xc, Yc] = [Hc, Vc] = meshgrid(Bc, Bc);
	Hxc = Vyc = Hc.*0 + 1;
	Hxyc = Vxyc = Hyc = Vxc = Hc.*0;                
	t
	end;
	
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




if(mod(ctr, 100)==0)
        
        t;

		%using sampling grid
		__x = H2dw(Xf, Yf, Hf, Hxf, Hyf, Hxyf, Xs, Ys);
		__y = H2dw(Xf, Yf, Vf, Vxf, Vyf, Vxyf, Xs, Ys);
        
end


ctr = ctr+1;

end;

