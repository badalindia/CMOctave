N = 64;
xDelta = 1;
x = linspace(0, xDelta, N+1);
h = x(2) - x(1);
[X, Y] = meshgrid(x);

%fourier frequencies
K = [0:N/2, -N/2+1:-1];
[kx, ky] = meshgrid(K, K);

%initial vorticity
W = sin(2*pi*X).*sin(2*pi*Y);

%obstacle function
[Pc00, Pc01, Pc10, Pc11] = obs_funct2(X, Y, h);

%constants
dt = 10^-3;
eta = 10^-2;
eh = 10^-6;

t=0;
tf = 1000;
while t<tf


	%Calculating velocity from vorticity
		W_ = fft2(W(1:N, 1:N));
			
		psi_ = -W_./((kx.^2 + ky.^2)*(2*pi/xDelta)^2);
		psi_(1,1) = 0;
		psix_ = (2*pi/xDelta)*1i.*psi_.*kx;
		psiy_ = (2*pi/xDelta)*1i.*psi_.*ky;
		psixy_ = -((2*pi/xDelta).^2)*psi_.*kx.*ky;

		psi = real( ifft2(psi_) )([1:N, 1], [1:N, 1]);
		psix = real( ifft2(psix_) )([1:N, 1], [1:N, 1]);
		psiy = real( ifft2(psiy_) )([1:N, 1], [1:N, 1]);
		psixy = real( ifft2(psixy_) )([1:N, 1], [1:N, 1]);

		%velocity on fine grid
		u = -real( ifft2(psiy_) )([1:N, 1], [1:N, 1]);
		v =  real( ifft2(psix_) )([1:N, 1], [1:N, 1]);


	%advection using fpi
	[x0_xpyp, y0_xpyp] = fpi(X, Y, X+eh, Y+eh, t, dt, psi, psix, psiy, psixy);
	[x0_xpyn, y0_xpyn] = fpi(X, Y, X+eh, Y-eh, t, dt, psi, psix, psiy, psixy);
	[x0_xnyp, y0_xnyp] = fpi(X, Y, X-eh, Y+eh, t, dt, psi, psix, psiy, psixy);
	[x0_xnyn, y0_xnyn] = fpi(X, Y, X-eh, Y-eh, t, dt, psi, psix, psiy, psixy);
	
	
	H0_xpyp  = H2dw(X, Y, H, Hx, Hy, Hxy, x0_xpyp, y0_xpyp);
	H0_xpyn  = H2dw(X, Y, H, Hx, Hy, Hxy, x0_xpyn, y0_xpyn);
	H0_xnyp  = H2dw(X, Y, H, Hx, Hy, Hxy, x0_xnyp, y0_xnyp);
	H0_xnyn  = H2dw(X, Y, H, Hx, Hy, Hxy, x0_xnyn, y0_xnyn);
	
	V0_xpyp  = H2dw(X, Y, V, Vx, Vy, Vxy, x0_xpyp, y0_xpyp);
	V0_xpyn  = H2dw(X, Y, V, Vx, Vy, Vxy, x0_xpyn, y0_xpyn);
	V0_xnyp  = H2dw(X, Y, V, Vx, Vy, Vxy, x0_xnyp, y0_xnyp);
	V0_xnyn  = H2dw(X, Y, V, Vx, Vy, Vxy, x0_xnyn, y0_xnyn);
	
	H = (H0_xpyp + H0_xnyn + H0_xpyn + H0_xnyp)/4;
	Hx = (H0_xpyp - H0_xnyn + H0_xpyn - H0_xnyp)/(4*eh);
	Hy = (H0_xpyp - H0_xnyn - H0_xpyn + H0_xnyp)/(4*eh);
	Hxy = (H0_xpyp + H0_xnyn - H0_xpyn - H0_xnyp)/(4*eh*eh);
	
	V = (V0_xpyp + V0_xnyn + V0_xpyn + V0_xnyp)/4;
	Vx = (V0_xpyp - V0_xnyn + V0_xpyn - V0_xnyp)/(4*eh);
	Vy = (V0_xpyp - V0_xnyn - V0_xpyn + V0_xnyp)/(4*eh);
	Vxy = (V0_xpyp + V0_xnyn - V0_xpyn - V0_xnyp)/(4*eh*eh);


	%adding penalty function
	du = - dt/eta*(Pc00.*u + Pc01.*v);
	dv = - dt/eta*(Pc10.*u + Pc11.*v);
	
	du_ = fft2(du([1:N], [1:N]));
	dv_ = fft2(dv([1:N], [1:N]));
	
	du_dy = (2*pi/xDelta)*1i*du_.*ky;
	dv_dx = (2*pi/xDelta)*1i*dv_.*kx;
	
	dW_ = dv_dx - du_dy;
	
	dW = real( ifft2(dW_)  )([1:N, 1],[1:N, 1]);
	
	W = W + dW;
	
	t = t+dt;
end;


th = linspace(0,2*pi,100)'; 
circsx = 0.15.*cos(th) + 0.5; 
circsy = 0.15.*sin(th) + 0.5; 

quiver(X, Y, u, v);
axis([0, 1, 0, 1], 'square');
	hold on;
	plot(circsx, circsy);
	drawnow;
	hold off;
print('staticV.png', '-r200');
