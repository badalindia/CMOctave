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
dL = 1;
[x, y] = meshgrid(linspace(0, dL, Nf+1), linspace(0, dL, Nf+1));
[Iy Ix] = size(x);
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
[Xs, Ys] = meshgrid(linspace(0, dL, Np+1));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Fw=@(x, y) sin(2*pi*x).*sin(2*pi*y);

%RHO		= 	@(x, y) exp(-12000*(((x-0.53).^2 + (y-0.45).^2).^2));
%T  		= 	@(x, y) exp(-100*((x-0.47).^2 + (y-0.3).^2));
coeff = 1.5;
Tamb = 0;

alpha = 0;
beta  = 0.4;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Initial parameters
W0 = Wx0 = Wy0 = Wxy0 = zeros(Nf+1, Nf+1);
t0 = 0;
tf = 2*pi;
dt = 2*10^-3;
eh = 10^-4;	
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


%forward particles
Xfp = Xf;
Yfp = Yf;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%some control parameters
ifPlot = true;
nC = 64;		%Number of contour bands to be displayed
gridShiftN = 2;	
plotctr = 1;	%No of step to wait before each plotting
t = t0;			
i = j = [2:Nc];
ctr = 0;
Res = '-r100';
target = "Smoke";

G = W0;
%G = 2*pi*cos(t).*(cos(2*pi*x).*sin(pi*y).^2 + cos(2*pi*y).*sin(pi*x).^2);

if(ifPlot && mod(ctr, plotctr)==0)

	myplot(2, Xs, Ys, 1-RHO(Xs, Ys, coeff), nC, gray, [0, 1], 
	strcat('Images/', target, '/rho/RHO_',num2str(Nc), '_',num2str(Nf), '_',num2str(Np), '_', num2str(ctr),'.png'), Res);


	myplot(3, Xs, Ys, T(Xs, Ys, coeff), nC, hot, [0, 1], 
	strcat('Images/', target, '/t/Temp_',num2str(Nc), '_',num2str(Nf), '_',num2str(Np), '_', num2str(ctr),'.png'), Res);
	
	
	myplot(4, Xf, Yf, G, nC, jet, 0, 
	strcat('Images/', target, '/w/Vorticity_',num2str(Nc), '_',num2str(Nf), '_',num2str(Np), '_', num2str(ctr),'.png'), Res);
	
end;

ctr = ctr+1;

while t<tf;
		
		%evaluating velocity field from F
		
		%Temporary vorticity field
			%GTemp = 2*pi*cos(t).*(cos(2*pi*x).*sin(pi*y).^2 + cos(2*pi*y).*sin(pi*x).^2);
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
        
        Ac = real(ifft2(A_c))([1:Nc, 1], [1:Nc, 1]);
		Axc = real(ifft2(Ax_c))([1:Nc, 1], [1:Nc, 1]);
		Ayc = real(ifft2(Ay_c))([1:Nc, 1], [1:Nc, 1]);
		Axyc = real(ifft2(Axy_c))([1:Nc, 1], [1:Nc, 1]);
		
		%velocity on fine grid
		Uxf = -real( ifft2(psify_) )([1:Nf, 1], [1:Nf, 1]);
		Uyf =  real( ifft2(psifx_) )([1:Nf, 1], [1:Nf, 1]);
		
	dt = min(0.01, 0.5*h/max(max(sqrt(Uxf.^2 + Uyf.^2))) );
		
	%advection using fpi
	[x0_xpyp, y0_xpyp] = fpi(Xc, Yc, Xc+eh, Yc+eh, t, dt, Ac, Axc, Ayc, Axyc);
	[x0_xpyn, y0_xpyn] = fpi(Xc, Yc, Xc+eh, Yc-eh, t, dt, Ac, Axc, Ayc, Axyc);
	[x0_xnyp, y0_xnyp] = fpi(Xc, Yc, Xc-eh, Yc+eh, t, dt, Ac, Axc, Ayc, Axyc);
	[x0_xnyn, y0_xnyn] = fpi(Xc, Yc, Xc-eh, Yc-eh, t, dt, Ac, Axc, Ayc, Axyc);
	
	
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
	
	%updating forward particles
	[psix_fp, psiy_fp] = H2dw_dxdy(Xc, Yc, Ac, Axc, Ayc, Axyc, Xfp, Yfp);
	Uxfp = -psiy_fp;
	Uyfp =  psix_fp;
	
	Xfp = Xfp + Uxfp*dt;
	Yfp = Yfp + Uyfp*dt;
	
	%temporarily composed map
	[Hf_t, Hxf_t, Hyf_t, Hxyf_t, Vf_t, Vxf_t, Vyf_t, Vxyf_t] = composeMaps(	Xf, Yf, xDelta, yDelta,
																Xc, Yc, H, Hx, Hy, Hxy, V, Vx, Vy, Vxy, 
																Xf, Yf, Hf, Hxf, Hyf, Hxyf, Vf, Vxf, Vyf, Vxyf);
		
		
	_x = H2dw(Xf, Yf, Hf_t, Hxf_t, Hyf_t, Hxyf_t, x, y);
	_y = H2dw(Xf, Yf, Vf_t, Vxf_t, Vyf_t, Vxyf_t, x, y);
	
	
		%updating source term
		_x_xn = H2dw(Xf, Yf, Hf_t, Hxf_t, Hyf_t, Hxyf_t, Xfp-eh, Yfp);
		_y_xn = H2dw(Xf, Yf, Vf_t, Vxf_t, Vyf_t, Vxyf_t, Xfp-eh, Yfp);
		_x_xp = H2dw(Xf, Yf, Hf_t, Hxf_t, Hyf_t, Hxyf_t, Xfp+eh, Yfp);
		_y_xp = H2dw(Xf, Yf, Vf_t, Vxf_t, Vyf_t, Vxyf_t, Xfp+eh, Yfp);
		
		dW = (-alpha*(RHO(_x_xp, _y_xp, coeff) - RHO(_x_xn, _y_xn, coeff))./(2*eh) + beta*(T(_x_xp, _y_xp, coeff) - T(_x_xn, _y_xn, coeff))./(2*eh))*dt;
		
		
		#{
		f_xp = beta*(T(_x_xp, _y_xp, coeff) - Tamb)./(RHO(_x_xp, _y_xp, coeff) + 10^-5);
		f_xn = beta*(T(_x_xn, _y_xn, coeff) - Tamb)./(RHO(_x_xn, _y_xn, coeff) + 10^-5);
		dW = dt.*(f_xp - f_xn)./(2*eh);
		#}
		
		W0 = W0 + dW;
		
		W0_   = fft2(W0([1:Nf], [1:Nf]));
		Wx0_  =  (2*pi/yDelta)*1i*W0_.*kx;
		Wy0_  =  (2*pi/yDelta)*1i*W0_.*ky;
		Wxy0_ = -(2*pi/yDelta)^2 *W0_.*kx.*ky;
		
		Wx0  = real( ifft2(Wx0_)  )([1:Nf, 1],[1:Nf, 1]);
		Wy0  = real( ifft2(Wy0_)  )([1:Nf, 1],[1:Nf, 1]);
		Wxy0 = real( ifft2(Wxy0_) )([1:Nf, 1],[1:Nf, 1]);	
		
		G = H2dw(Xf, Yf, W0, Wx0, Wy0, Wxy0, _x, _y);
			
	t = t+dt;
	
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
	end;
	
	if(ifPlot && mod(ctr, plotctr)==0)
	
		dt
	
			%using sampling grid
			__x = H2dw(Xf, Yf, Hf_t, Hxf_t, Hyf_t, Hxyf_t, Xs, Ys);
			__y = H2dw(Xf, Yf, Vf_t, Vxf_t, Vyf_t, Vxyf_t, Xs, Ys);
		
		G1 = 1-RHO(__x, __y, coeff);
		G2 = T(__x, __y, coeff);
		G3 = H2dw(Xf, Yf, W0, Wx0, Wy0, Wxy0, __x, __y);	
		%G3 = bilinearInterpolation(x, y, __x, __y, W0);
				
		myplot(2, Xs, Ys, G1, nC, gray, [0, 1], 
		strcat('Images/', target, '/rho/RHO_',num2str(Nc), '_',num2str(Nf), '_',num2str(Np), '_', num2str(ctr),'.png'), Res);


		myplot(3, Xs, Ys, G2, nC, hot, [0, 1], 
		strcat('Images/', target, '/t/Temp_',num2str(Nc), '_',num2str(Nf), '_',num2str(Np), '_', num2str(ctr),'.png'), Res);
		
		
		myplot(4, Xs, Ys, G3, nC, jet, 0, 
		strcat('Images/', target, '/w/Vorticity_',num2str(Nc), '_',num2str(Nf), '_',num2str(Np), '_', num2str(ctr),'.png'), Res);
		
			
	end;
	
	if(mod(ctr, 10)==0)
		t
	end
	
	ctr = ctr+1;
	
end
