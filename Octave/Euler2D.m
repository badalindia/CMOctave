clear;

page_output_immediately(true);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Initial conditions
Nc = 32;
Nf = 256;
Np = 512;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%initialization of reference grid
dL = 2*pi;
[x, y] = meshgrid(linspace(0, dL, Nf+1), linspace(0, dL, Nf+1));
[Iy Ix] = size(x)
xMin = x(1, 1);
xMax = x(1, Ix);
xDelta = xMax - xMin;
yMin = y(1, 1);
yMax = y(Iy, 1);
yDelta = yMax - yMin;
h = x(1, 2) - x(1, 1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sampling grid
[Xs, Ys] = meshgrid(linspace(0, dL, Np));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Initial vorticity function
Fw = @(x, y) -sin(x*pi/dL).^2.*sin(y*pi/dL).^2.*(exp(-((x - 0.5*dL).^2 + (y - dL/3).^2)*5) + exp(-((x - 0.5*dL).^2 + (y - 2*dL/3).^2)*5));
%Fw = @(x, y) sin(x).*sin(y);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Initial parameters
F = Fw(x, y);
peak = max(max(abs(F)))
t0 = 0;
tf = 2000;
ifPlot = true;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%diffeomorphism grids

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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%extending (warping) domain on corners for periodic conditions
[IXc_w, IYc_w] = meshgrid([1:Nc+3], [1:Nc+3]);	%w for warped

ifXcWarped = ((IXc_w-1)>Nc+1) - ((IXc_w-1)<1);
ifYcWarped = ((IYc_w-1)>Nc+1) - ((IYc_w-1)<1);

XcMax = Xc(1, Nc+1);
XcMin = Xc(1, 1);
XcDelta = XcMax - XcMin;
YcMax = Yc(Nc+1, 1);
YcMin = Yc(1, 1);
YcDelta = YcMax - YcMin;

IXc_w = (IXc_w-1) - (Nc)*((IXc_w-1)>Nc+1) +  (Nc)*((IXc_w-1)<1);
IYc_w = (IYc_w-1) - (Nc)*((IYc_w-1)>Nc+1) +  (Nc)*((IYc_w-1)<1);

Ic_w = (IXc_w-1).*(Nc+1) + IYc_w;
Xc_w = Xc(Ic_w) + ifXcWarped.*XcDelta;
Yc_w = Yc(Ic_w) + ifYcWarped.*YcDelta;

[IXf_w, IYf_w] = meshgrid([1:Nf+3], [1:Nf+3]);	%w for warped

ifXfWarped = ((IXf_w-1)>Nf+1) - ((IXf_w-1)<1);
ifYfWarped = ((IYf_w-1)>Nf+1) - ((IYf_w-1)<1);

XfMax = Xf(1, Nf+1);
XfMin = Xf(1, 1);
XfDelta = XfMax - XfMin;
YfMax = Yf(Nf+1, 1);
YfMin = Yf(1, 1);
YfDelta = YfMax - YfMin;

IXf_w = (IXf_w-1) - (Nf)*((IXf_w-1)>Nf+1) +  (Nf)*((IXf_w-1)<1);
IYf_w = (IYf_w-1) - (Nf)*((IYf_w-1)>Nf+1) +  (Nf)*((IYf_w-1)<1);

If_w = (IXf_w-1).*(Nf+1) + IYf_w;
Xf_w = Xf(If_w) + ifXfWarped.*XfDelta;
Yf_w = Yf(If_w) + ifYfWarped.*YfDelta;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%some control parameters
nC = 100;		%Number of contour bands to be displayed
gridShiftN = 5;	
plotctr = 50;	%No of step to wait before each plotting
t = t0;			
i = j = [2:Nc];
ctr = 0;

G = F;
	
if(ifPlot && mod(ctr, plotctr)==0)

		%using sampling grid
		Gtemp = Fw(Xs, Ys);
		
	
	contourf(Xs, Ys, Gtemp, nC);
	shading flat;
	axis([xMin, xMax, yMin, yMax], "square");
	caxis([-peak, 0]);
	shading flat;	
	drawnow;
	print( strcat('Images/Eule2D_',num2str(Nc), '_',num2str(Nf), '_', num2str(ctr),'.png'), '-r200');
end;

ctr = ctr+1;
dt = hc/2

while t<tf;
		
	if(t+dt>tf)
		dt = tf - t;
	end;	
		
		%evaluating velocity field from F
		
		%Temporary vorticity field
		F_ = fft2(G(1:Nf, 1:Nf));
        
        K = [0:Nf/2, -Nf/2+1:-1];
		[kx, ky] = meshgrid(K, K);
        
        psif_ = -F_./((kx.^2 + ky.^2)*(2*pi/xDelta)^2);
        psif_(1,1) = 0;
        psifx_ = (2*pi/xDelta)*1i.*psif_.*kx;
        psify_ = (2*pi/xDelta)*1i.*psif_.*ky;
        psifxy_ = -((2*pi/xDelta).^2)*psif_.*kx.*ky;
        
        %mask
        kf2c = [1:Nc/2+1, Nf-(Nc/2-1)+1:Nf];
        M = zeros(Nc, Nc) + (Nc/Nf)^2;
        M(Nc/2+1, :) = M(Nc/2+1, :)*2;
        M(:, Nc/2+1) = M(:, Nc/2+1)*2;
        
        if(Nc==Nf)
            printf('Mask is not designed for Nc=Nf');
        end
        
        %on coarser grid
        
        A_c = psif_(kf2c, kf2c).*M;
        Ax_c = psifx_(kf2c, kf2c).*M;
        Ay_c = psify_(kf2c, kf2c).*M;
        Axy_c = psifxy_(kf2c, kf2c).*M;
        
        %Kc = [0:Nc/2, -Nc/2+1:-1];
        %[kxc, kyc] = meshgrid(Kc, Kc);
        	
        Ac = real(ifft2(A_c))([1:Nc, 1], [1:Nc, 1]);
		Axc = real(ifft2(Ax_c))([1:Nc, 1], [1:Nc, 1]);
		Ayc = real(ifft2(Ay_c))([1:Nc, 1], [1:Nc, 1]);
		Axyc = real(ifft2(Axy_c))([1:Nc, 1], [1:Nc, 1]);
        
		Ac_w = Ac(Ic_w);
		Axc_w = Axc(Ic_w);
		Ayc_w = Ayc(Ic_w);
		Axyc_w = Axyc(Ic_w);
		
	%backward integration
	eh = 10^-4;
	
    %advection using Euler
	[x0_xpyp, y0_xpyp] = Advection_Euler1D(Xc_w, Yc_w, Xc+eh, Yc+eh, dt, Ac_w, Axc_w, Ayc_w, Axyc_w);
	[x0_xpyn, y0_xpyn] = Advection_Euler1D(Xc_w, Yc_w, Xc+eh, Yc-eh, dt, Ac_w, Axc_w, Ayc_w, Axyc_w);
	[x0_xnyp, y0_xnyp] = Advection_Euler1D(Xc_w, Yc_w, Xc-eh, Yc+eh, dt, Ac_w, Axc_w, Ayc_w, Axyc_w);
	[x0_xnyn, y0_xnyn] = Advection_Euler1D(Xc_w, Yc_w, Xc-eh, Yc-eh, dt, Ac_w, Axc_w, Ayc_w, Axyc_w);
	%//////////////////////////////
    
    #{
	%advection using fpi
	[x0_xpyp, y0_xpyp] = fpi(Xc_w, Yc_w, Xc+eh, Yc+eh, t, dt, Ac_w, Axc_w, Ayc_w, Axyc_w);
	[x0_xpyn, y0_xpyn] = fpi(Xc_w, Yc_w, Xc+eh, Yc-eh, t, dt, Ac_w, Axc_w, Ayc_w, Axyc_w);
	[x0_xnyp, y0_xnyp] = fpi(Xc_w, Yc_w, Xc-eh, Yc+eh, t, dt, Ac_w, Axc_w, Ayc_w, Axyc_w);
	[x0_xnyn, y0_xnyn] = fpi(Xc_w, Yc_w, Xc-eh, Yc-eh, t, dt, Ac_w, Axc_w, Ayc_w, Axyc_w);
	%//////////////////////////////
	#}
	
	H0_xpyp  = H2d_Warped(Xc, Yc, H, Hx, Hy, Hxy, xDelta, 0, x0_xpyp, y0_xpyp);
	H0_xpyn  = H2d_Warped(Xc, Yc, H, Hx, Hy, Hxy, xDelta, 0, x0_xpyn, y0_xpyn);
	H0_xnyp  = H2d_Warped(Xc, Yc, H, Hx, Hy, Hxy, xDelta, 0, x0_xnyp, y0_xnyp);
	H0_xnyn  = H2d_Warped(Xc, Yc, H, Hx, Hy, Hxy, xDelta, 0, x0_xnyn, y0_xnyn);
	
	V0_xpyp  = H2d_Warped(Xc, Yc, V, Vx, Vy, Vxy, 0, yDelta, x0_xpyp, y0_xpyp);
	V0_xpyn  = H2d_Warped(Xc, Yc, V, Vx, Vy, Vxy, 0, yDelta, x0_xpyn, y0_xpyn);
	V0_xnyp  = H2d_Warped(Xc, Yc, V, Vx, Vy, Vxy, 0, yDelta, x0_xnyp, y0_xnyp);
	V0_xnyn  = H2d_Warped(Xc, Yc, V, Vx, Vy, Vxy, 0, yDelta, x0_xnyn, y0_xnyn);
	
	H = (H0_xpyp + H0_xnyn + H0_xpyn + H0_xnyp)/4;
	Hx = (H0_xpyp - H0_xnyn + H0_xpyn - H0_xnyp)/(4*eh);
	Hy = (H0_xpyp - H0_xnyn - H0_xpyn + H0_xnyp)/(4*eh);
	Hxy = (H0_xpyp + H0_xnyn - H0_xpyn - H0_xnyp)/(4*eh*eh);
		
	V = (V0_xpyp + V0_xnyn + V0_xpyn + V0_xnyp)/4;
	Vx = (V0_xpyp - V0_xnyn + V0_xpyn - V0_xnyp)/(4*eh);
	Vy = (V0_xpyp - V0_xnyn - V0_xpyn + V0_xnyp)/(4*eh);
	Vxy = (V0_xpyp + V0_xnyn - V0_xpyn - V0_xnyp)/(4*eh*eh);
	
	%temporarily composed map
	[Hf_t, Hxf_t, Hyf_t, Hxyf_t, Vf_t, Vxf_t, Vyf_t, Vxyf_t] = composeMaps(	Xf, Yf, xDelta, yDelta,
																Xc, Yc, H, Hx, Hy, Hxy, V, Vx, Vy, Vxy, 
																Xf, Yf, Hf, Hxf, Hyf, Hxyf, Vf, Vxf, Vyf, Vxyf);
																
	%advected field
	_x = H2d_Warped(Xf, Yf, Hf_t, Hxf_t, Hyf_t, Hxyf_t, xDelta, 0, x, y);
	_y = H2d_Warped(Xf, Yf, Vf_t, Vxf_t, Vyf_t, Vxyf_t, 0, yDelta, x, y);
	
	G = Fw(_x, _y);
	
	if(ifPlot && mod(ctr, plotctr)==0)
	
			%using sampling grid
			__x = H2d_Warped(Xf, Yf, Hf_t, Hxf_t, Hyf_t, Hxyf_t, xDelta, 0, Xs, Ys);
			__y = H2d_Warped(Xf, Yf, Vf_t, Vxf_t, Vyf_t, Vxyf_t, 0, yDelta, Xs, Ys);
			Gtemp = Fw(__x, __y);
			
	
		contourf(Xs, Ys, Gtemp, nC);
		axis([xMin, xMax, yMin, yMax], "square");
		caxis([-peak, 0]);
		shading flat;
		drawnow;
		print( strcat('Images/Euler2D_',num2str(Nc), '_',num2str(Nf), '_', num2str(ctr),'.png'), '-r200');
		shading flat;
		t
	end;
	
		%grid composition and resetting
	if(mod(ctr,gridShiftN)==0)
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
	
	t = t+dt;
	ctr = ctr+1;
end

