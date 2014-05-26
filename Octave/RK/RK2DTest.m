%tests global accuracy of RK method using exact solutions
clear all;

%test function definitoin
	%dx/dt = u(x(t), t)
u = @(x, y, t) -y + 2*t.*(cos(t) - x);
v = @(x, y, t) x + exp(-t.*t);

xExact = @(t) cos(t) - exp(-t.*t);
yExact = @(t) sin(t);

%test parameters
t0 = 0;			%initial time
tf = pi/2;		%final time
n = [2:10];		%time domain will divide into 2^n equal step sizes


x0 = xExact(t0);
y0 = yExact(t0);
xfExact = xExact(tf);
yfExact = yExact(tf);
N = 2.^n;
Dt = (tf - t0)./N;
I = length(n);

e = zeros(I, 1);	%global error

for i=[1:I]

	t=0;
	x = x0;
	y = y0;
	dt = Dt(i);
	
	while t<tf - 0.5*dt;	
		[x, y] = RK3_2D(x, y, t, dt, u, v);
		t = t+dt;
	end	
	
	e(i) =  norm([xfExact-x, yfExact-y], 2);
end;

plot(log(Dt), log(e), '-o');
Accuracy = polyfit(log(Dt'), log(e), 1)(1)

