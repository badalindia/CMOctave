clear;

page_output_immediately(true);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Initial conditions
Nc = 64;
Nf = 256;
Np = 512;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Grid
dL = 1;
x = linspace(0, dL, Nc+1);
h = x(2) - x(1);
[X, Y] = meshgrid(x);

K = [0:Nc/2, -Nc/2+1:-1];
[kx, ky] = meshgrid(K, K);
ksquare = kx.^2 + ky.^2;

%mask
kf2c = [1:Nc/2+1, Nf-(Nc/2-1)+1:Nf];
M = zeros(Nc, Nc) + (Nc/Nf)^2;
M(Nc/2+1, :) = M(Nc/2+1, :)*2;
M(:, Nc/2+1) = M(:, Nc/2+1)*2;


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



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Initial fields
%W = exp(-((X-dL/2).^2+(Y-dL/2+dL/8).^2)/(0.01))+exp(-((X-dL/2).^2+(Y-dL/2-dL/8).^2)/(0.01))-0.5*exp(-((X-dL/2-dL/8).^2+(Y-dL/2-dL/8).^2)/(0.02));
W = X.*0;

        %calculation of initial velocity
        W_ = fft2(W(1:Nc, 1:Nc));

        psi_ = -W_./((kx.^2 + ky.^2)*(2*pi/dL)^2);
        psi_(1,1) = 0;
        psix_ = (2*pi/dL)*1i.*psi_.*kx;
        psiy_ = (2*pi/dL)*1i.*psi_.*ky;
        psixy_ = -((2*pi/dL).^2)*psi_.*kx.*ky;
        
        U_ = -psiy_;
        V_ =  psix_;        

%penalty
P = obs_functNS(X, Y, h);


nu = 10^-4;  % viscosity        
coeff = 1.5;
alpha = 0;
beta  = 0.4;
eh = 10^-6;
gridShiftN = 5;
eta = 10^-2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plotting parameters
plotCtr = 200;
target = "NavierStokes_Smoke_Obs";
Res = '-r150';
nC = 50;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


t = 0;
dt = 5*10^-4;
tf = 40;
ctr = 1;
ctr2 = 0;

t1 = clock();
while t<tf;
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
	

	 %calculation of buoyancy force
	_x0 = H2dw(Xc, Yc, Hc, Hxc, Hyc, Hxyc, Xf, Yf);
	_y0 = H2dw(Xc, Yc, Vc, Vxc, Vyc, Vxyc, Xf, Yf);
	
	_x = H2dw(Xf, Yf, Hf, Hxf, Hyf, Hxyf, _x0, _y0);
	_y = H2dw(Xf, Yf, Vf, Vxf, Vyf, Vxyf, _x0, _y0);
	
	Fb  = -alpha*RHO(_x, _y, coeff) + beta*T(_x, _y, coeff);
	Fb_ = fft2(Fb([1:Nf], [1:Nf]));
	Fbc_ = Fb_(kf2c, kf2c).*M;	
	
	Fox = - P.*U/eta;
	Foy = - P.*V/eta;
	
	Fox_ = fft2(Fox([1:Nc], [1:Nc]));
	Foy_ = fft2(Foy([1:Nc], [1:Nc]));
		
	%Cranck Nicolson time marching
	Ustar_ = (U_.*(1 - nu*dt/2*(2*pi/dL)^2*ksquare) + dt*(-ConvX_        + Fox_ ))./(1 + nu*dt/2*(2*pi/dL)^2*ksquare);
	Vstar_ = (V_.*(1 - nu*dt/2*(2*pi/dL)^2*ksquare) + dt*(-ConvY_ + Fbc_ + Foy_ ))./(1 + nu*dt/2*(2*pi/dL)^2*ksquare);
	
%step 2, calculating pressure (solving Poisson's equation)

	divU_ = (2*pi/dL)*1i.*(Ustar_.*kx + Vstar_.*ky);
	P_ = -divU_./((2*pi/dL)^2*dt*ksquare);
	P_(1) = 0;
	Px_ = (2*pi/dL)*1i.*P_.*kx;
	Py_ = (2*pi/dL)*1i.*P_.*ky;
	
%step 3, calculating Unew
	
	U_ = Ustar_ - dt*Px_;
	V_ = Vstar_ - dt*Py_;
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%advection of chi
	U0_   =  U_;
	U0x_  =  (2*pi/dL)*1i.*U_.*kx;
	U0y_  =  (2*pi/dL)*1i.*U_.*ky;
	U0xy_ = -(2*pi/dL).^2.*U_.*kx.*ky;
	
	V0_   =  V_;
	V0x_  =  (2*pi/dL)*1i.*V_.*kx;
	V0y_  =  (2*pi/dL)*1i.*V_.*ky;
	V0xy_ = -(2*pi/dL).^2.*V_.*kx.*ky;
	
	U0   = real( ifft2(U0_  ) )([1:Nc, 1], [1:Nc, 1]);
	U0x  = real( ifft2(U0x_ ) )([1:Nc, 1], [1:Nc, 1]);
	U0y  = real( ifft2(U0y_ ) )([1:Nc, 1], [1:Nc, 1]);
	U0xy = real( ifft2(U0xy_) )([1:Nc, 1], [1:Nc, 1]);
	
	V0   = real( ifft2(V0_  ) )([1:Nc, 1], [1:Nc, 1]);
	V0x  = real( ifft2(V0x_ ) )([1:Nc, 1], [1:Nc, 1]);
	V0y  = real( ifft2(V0y_ ) )([1:Nc, 1], [1:Nc, 1]);
	V0xy = real( ifft2(V0xy_) )([1:Nc, 1], [1:Nc, 1]);
	
	%advection using Euler
	[x0_xpyp, y0_xpyp] = AdvectionUV_Euler1D(Xc, Yc, Xc+eh, Yc+eh, dt, U0, U0x, U0y, U0xy, V0, V0x, V0y, V0xy);
	[x0_xpyn, y0_xpyn] = AdvectionUV_Euler1D(Xc, Yc, Xc+eh, Yc-eh, dt, U0, U0x, U0y, U0xy, V0, V0x, V0y, V0xy);
	[x0_xnyp, y0_xnyp] = AdvectionUV_Euler1D(Xc, Yc, Xc-eh, Yc+eh, dt, U0, U0x, U0y, U0xy, V0, V0x, V0y, V0xy);
	[x0_xnyn, y0_xnyn] = AdvectionUV_Euler1D(Xc, Yc, Xc-eh, Yc-eh, dt, U0, U0x, U0y, U0xy, V0, V0x, V0y, V0xy);
	
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
	
	
	%grid composition and resetting
	if(mod(ctr,gridShiftN)==0)
	[Hf, Hxf, Hyf, Hxyf, Vf, Vxf, Vyf, Vxyf] = composeMaps(Xf, Yf, dL, dL,
												Xc, Yc, Hc, Hxc, Hyc, Hxyc, Vc, Vxc, Vyc, Vxyc, 
												Xf, Yf, Hf, Hxf, Hyf, Hxyf, Vf, Vxf, Vyf, Vxyf);

	[Xc, Yc] = [Hc, Vc] = meshgrid(Bc, Bc);
	Hxc = Vyc = Hc.*0 + 1;
	Hxyc = Vxyc = Hyc = Vxc = Hc.*0;                
	end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


t = t+dt;


if(mod(ctr, plotCtr)==0)
        
        t2 = clock();
        time = etime(t2, t1)/plotCtr
        dt
        t

	            %using sampling grid
                __x = H2dw(Xf, Yf, Hf, Hxf, Hyf, Hxyf, Xs, Ys);
                __y = H2dw(Xf, Yf, Vf, Vxf, Vyf, Vxyf, Xs, Ys);
                
                G1 = 1-RHO(__x, __y, coeff);
                G2 = T(__x, __y, coeff);
        
        myplot2(1, Xs, Ys, G1, nC, gray, [0, 1], 
        strcat('Images/', target, '/rho/RHO_',num2str(Nc), '_',num2str(Nf), '_',num2str(Np), '_', num2str(ctr2),'.png'), Res);

		myplot2(2, Xs, Ys, G2, nC, hot, [0, 1], 
        strcat('Images/', target, '/t/Temp_',num2str(Nc), '_',num2str(Nf), '_',num2str(Np), '_', num2str(ctr2),'.png'), Res);
        
        ctr2 = ctr2+1;
        
        t1 = clock();
end

ctr = ctr+1;

end;

