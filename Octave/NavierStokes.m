%Navier Stokes solved using projection method

clear;

page_output_immediately(true);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Initial conditions
Nc = 32;
Nf = 64;
Np = 512;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Grid
dL = 1;
x = linspace(0, dL, Nc+1);
[X, Y] = meshgrid(x);

K = [0:Nc/2, -Nc/2+1:-1];
[kx, ky] = meshgrid(K, K);
ksquare = kx.^2 + ky.^2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Initial fields
%W = exp(-((X-pi).^2+(Y-pi+pi/4).^2)/(0.2))+exp(-((X-pi).^2+(Y-pi-pi/4).^2)/(0.2))-0.5*exp(-((X-pi-pi/4).^2+(Y-pi-pi/4).^2)/(0.4));
W = exp(-((X-dL/2).^2+(Y-dL/2+dL/8).^2)/(0.01))+exp(-((X-dL/2).^2+(Y-dL/2-dL/8).^2)/(0.01))-0.5*exp(-((X-dL/2-dL/8).^2+(Y-dL/2-dL/8).^2)/(0.02));


	%calculation of initial velocity
	W_ = fft2(W(1:Nc, 1:Nc));

	psif_ = -W_./((kx.^2 + ky.^2)*(2*pi/dL)^2);
	psif_(1,1) = 0;
	psifx_ = (2*pi/dL)*1i.*psif_.*kx;
	psify_ = (2*pi/dL)*1i.*psif_.*ky;
	psifxy_ = -((2*pi/dL).^2)*psif_.*kx.*ky;
	
	U_ = -psify_;
	V_ =  psifx_;
	
		
		#{
		%temporary initialization
		U = sin(2*pi*X) + 2*sin(pi*X).^2.*cos(2*pi*Y);
		V = cos(2*pi*Y) + 2*cos(pi*Y).^2.*sin(2*pi*X);
		
		U_ = fft2(U(1:Nc, 1:Nc));
		V_ = fft2(V(1:Nc, 1:Nc));
		#}
		
			
	
nu=10^-4;  % viscosity	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plotting parameters
plotCtr = 100;
target = "NavierStokes";
Res = '-r150';
nC  = 50;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



t = 0;
dt = 10^-2;
tf = 2000*dt;
ctr = 1;
ctr2 = 0;


while t<tf-dt*0.1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% projection method


%step 1, calculating U*

	%calculating the convective flux (non-linear term)
	Ux_ = (2*pi/dL)*1i.*U_.*kx;
	Uy_ = (2*pi/dL)*1i.*U_.*ky;
	Vx_ = (2*pi/dL)*1i.*V_.*kx;
	Vy_ = (2*pi/dL)*1i.*V_.*ky;
	
	U = real( ifft2(U_) )([1:Nc, 1], [1:Nc, 1]);
	V = real( ifft2(V_) )([1:Nc, 1], [1:Nc, 1]);
	Ux = real( ifft2(Ux_) )([1:Nc, 1], [1:Nc, 1]);
	Uy = real( ifft2(Uy_) )([1:Nc, 1], [1:Nc, 1]);
	Vx = real( ifft2(Vx_) )([1:Nc, 1], [1:Nc, 1]);
	Vy = real( ifft2(Vy_) )([1:Nc, 1], [1:Nc, 1]);
	
	ConvX = (U.*Ux + V.*Uy);
	ConvY = (U.*Vx + V.*Vy);
	
	ConvX_ = fft2(ConvX([1:Nc], [1:Nc]));
	ConvY_ = fft2(ConvY([1:Nc], [1:Nc]));
	
	
		%Euler
		%Ustar_ = U_ + dt*(-ConvX_ - nu*ksquare.*U_);
		%Vstar_ = V_ + dt*(-ConvY_ - nu*ksquare.*V_);	
	
	%Cranck Nicolson	
	K = nu*dt/2*(2*pi/dL)^2*ksquare;
	
	Ustar_ = (U_.*(1 - K) + dt*(-ConvX_))./(1 + K);
	Vstar_ = (V_.*(1 - K) + dt*(-ConvY_))./(1 + K);
	
%step 2, calculating pressure (solving Poisson's equation)

	divU_ = (2*pi/dL)*1i.*(Ustar_.*kx + Vstar_.*ky);
	P_ = -divU_./((2*pi/dL)^2*dt*ksquare);
	P_(1) = 0;
	Px_ = (2*pi/dL)*1i.*P_.*kx; 
	Py_ = (2*pi/dL)*1i.*P_.*ky; 
	
	
%step 3, calculating Unew
	
	U_ = Ustar_ - dt*Px_;
	V_ = Vstar_ - dt*Py_;
	
t = t+dt;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(mod(ctr, plotCtr)==0)
	
	t

	Uy_ = (2*pi/dL)*1i.*U_.*ky;
	Vx_ = (2*pi/dL)*1i.*V_.*kx;
	
	W_ = Vx_ - Uy_;
	W = real( ifft2(W_) )([1:Nc, 1], [1:Nc, 1]); 
	
	myplot(1, X, Y, W, nC, jet, 0, 
	strcat('Images/', target, '/w/Vorticity_',num2str(Nc), '_',num2str(Nf), '_',num2str(Np), '_', num2str(ctr2),'.png'), Res);
	ctr2 = ctr2+1;
end

ctr = ctr+1;

end;
