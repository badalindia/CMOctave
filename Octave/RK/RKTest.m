%tests global accuracy of RK method using exact solutions
clear all;

%test function definitoin
	%dx/dt = u(x(t), t)
u = @(x, t) exp(-x)*sin(2*t);
xExact = @(t) log(1 + (sin(t)).^2);

%test parameters
t0 = 0;			%initial time
tf = pi/2;		%final time
n = [2:10];		%time domain will divide into 2^n equal step sizes


x0 = xExact(t0);
xfExact = xExact(tf);
N = 2.^n;
Dt = (tf - t0)./N;
I = length(n);

e = zeros(I, 1);	%global error

for i=[1:I]

	t=0;
	x = x0;
	dt = Dt(i);
	
	while t<tf - 0.5*dt;	
		x = RK3(x, t, dt, u);
		t = t+dt;
	end	
	
	e(i) =  abs(xfExact - x);	
end;

plot(log(Dt), log(e), '-o');
Accuracy = polyfit(log(Dt'), log(e), 1)(1)

