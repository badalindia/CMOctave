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
return;

%constants
dt = 10^-3;
eta = 10^-2;

t=0;
tf = 1;
while t<tf

	W_ = fft2(W(1:N, 1:N));
        
	psif_ = -W_./((kx.^2 + ky.^2)*(2*pi/xDelta)^2);
	psif_(1,1) = 0;
	psifx_ = (2*pi/xDelta)*1i.*psif_.*kx;
	psify_ = (2*pi/xDelta)*1i.*psif_.*ky;
	psifxy_ = -((2*pi/xDelta).^2)*psif_.*kx.*ky;

	%velocity on fine grid
	u = -real( ifft2(psify_) )([1:N, 1], [1:N, 1]);
	v =  real( ifft2(psifx_) )([1:N, 1], [1:N, 1]);

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
