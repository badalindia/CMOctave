function _Q = Diffeo1D(N, Xc, Q, U, t0, tf, dt, ifPlot);

%for fft operations
k = [0:N/2, -N/2+1:-1];
xDelta = max(Xc) - min(Xc);

%1D diffeo
Hc = Xc;
Hxc = Xc.*0 + 1;

t = t0;
ep = 10^(-6);


ctr = 0;


while t<tf  

	if(t+dt>tf)
		dt = tf - t;
	end;
	
	%calculating field at present point
	_x = H1dw(Xc, Hc, Hxc, Xc);
	Q1 = Q(_x);

	Q1_ = fft(Q1(1:N));
	U1_ = -Q1_./((k*2*pi/xDelta).^2);
	U1_(1) = 0;
	U1 = real(ifft(U1_))([1:N, 1]);
	U1 = -5*U1;
	
	U1_ = fft(U1(1:N));
	U1x_ = (1i*k*2*pi/xDelta).*U1_;
	U1x = real(ifft(U1x_))([1:N, 1]);
	
	Up = H1dw(Xc, U1, U1x, Xc+ep);
	Un = H1dw(Xc, U1, U1x, Xc-ep);
	
	x0p = Xc + ep - Up*dt;
	x0n = Xc - ep - Un*dt;
	
	Hp = H1dw(Xc, Hc, Hxc, x0p);
	Hn = H1dw(Xc, Hc, Hxc, x0n);
	
	Hc = (Hp+Hn)/2;
	Hxc = (Hp-Hn)/(2*ep);
	
	[Hc, Hxc] = refineMap(Xc, Hc, Hxc, mod(ctr, 2));
		
	t = t+dt;
	
	if ifPlot && mod(ctr, 10) == 0
	
	
		%_x = H1dw(Xc, Hc, Hxc, Xc);
		%_Q = Q(_x);
		%plot(Xc, _Q, '-or', Xc, Q(Xc), '-og');
		
		_x1 = linspace(0, 1, 255);
		_x2 = H1dw(Xc, Hc, Hxc, _x1);
		plot(_x1, _x2, '-r', Xc, Hc, 'ob');
		axis([0, 1, -0.5, 1.5]);
		%print(strcat("plots/plot_", num2str(ctr), ".png"), '-r100');
		
		drawnow;
	
	end
	ctr = ctr+1;
end;

ctr

_x = H1dw(Xc, Hc, Hxc, Xc);
_Q = Q(_x);
		
		

