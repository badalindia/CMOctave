function _Q = Diffeo1D(N, Xc, Q, U, t0, tf, dt, ifPlot);

%for fft operations
k = [0:N/2, -N/2+1:-1];
xDelta = max(Xc) - min(Xc);

%1D diffeo
Hc = Xc;
Hxc = Xc.*0 + 1;

t = t0;
ep = 10^(-6);

U1_p = Xc.*0;
U1x_p = Xc.*0;

while t<tf

	if(t+dt>tf)
		dt = tf - t;
	end;
	
		#{
		%advection using analytical velocity
		x0p = RK3(Xc+ep, t+dt, -dt, U);
		x0n = RK3(Xc-ep, t+dt, -dt, U);
		#}
		
			%calculating field at present point
			_x = H1dw(Xc, Hc, Hxc, Xc);
			Q1 = Q(_x);
		
			Q1_ = fft(Q1(1:N));
			U1_ = -Q1_./((k*2*pi/xDelta).^2);
			U1_(1) = 0;
			U1 = real(ifft(U1_))([1:N, 1]);
			U1 = 40*U1.*U(Xc, t);
			
			U1_ = fft(U1(1:N));
			U1x_ = (1i*k*2*pi/xDelta).*U1_;
			U1x = real(ifft(U1x_))([1:N, 1]);
			
			Up = H1dw(Xc, U1, U1x, Xc+ep);
			Un = H1dw(Xc, U1, U1x, Xc-ep);
			
				%simple euler moment
				#{
				x0p = Xc + ep - Up*dt;
				x0n = Xc - ep - Un*dt;
				#}
				
				%fixed point iteration
				x0p = fpi3U(Xc, Xc+ep, t, dt, U1, U1x, U1_p, U1x_p);
				x0n = fpi3U(Xc, Xc-ep, t, dt, U1, U1x, U1_p, U1x_p);
				U1_p = U1;
				U1x_p = U1x;
			
	
	Hp = H1dw(Xc, Hc, Hxc, x0p);
	Hn = H1dw(Xc, Hc, Hxc, x0n);
	
	Hc = (Hp+Hn)/2;
	Hxc = (Hp-Hn)/(2*ep);
		
	t = t+dt;
	
	if ifPlot
	
		_x = H1dw(Xc, Hc, Hxc, Xc);
		_Q = Q(_x);
		plot(Xc, _Q, '-or', Xc, Q(Xc), '-og');
		drawnow;
	end
	
end;

_x = H1dw(Xc, Hc, Hxc, Xc);
_Q = Q(_x);

