%tests global accuracy of RK method (without using exact solutions)
clear all;

%test function definitoin
	%dx/dt = u(x(t), t)
u = @(x, y, t) -y + 2*t.*(cos(t) - x);
v = @(x, y, t) x + exp(-t.*t) + y;

%test parameters
t0 = 0;			%initial time
tf = pi/2;		%final time
n = [5:14];		%time domain will divide into 2^n equal step sizes


x0 = 0;
y0 = 0;
N = 2.^n;
Dt = (tf - t0)./N;
I = length(n);

X = zeros(I, 1);	%global error
Y = zeros(I, 1);	%global error
e = zeros(I, 1);	%global error


for i=[1:I]

	t=0;
	x = x0;
	y = y0;
	dt = Dt(i);
	
	while t<tf - 0.5*dt;	
		[x, y] = RK2_2D(x, y, t, dt, u, v);
		t = t+dt;
	end	
	
	X(i) = x;
	Y(i) = y;
end;

e = ((X - X(I)).^2 + (Y - Y(I)).^2).^0.5;

plot(log(Dt), log(e), '-o');
Accuracy = polyfit(log(Dt(1:I-1)'), log(e(1:I-1)), 1)(1)

