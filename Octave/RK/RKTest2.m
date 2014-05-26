%tests global accuracy of RK method (without exact solutions)
clear all;

%test function definitoin
	%dx/dt = u(x(t), t)
u = @(x, t) exp(-x)*sin(2*t);

%test parameters
t0 = 0;			%initial time
tf = pi/2;		%final time
n = [2:10];		%time domain will divide into 2^n equal step sizes


x0 = 1;
N = 2.^n;
Dt = (tf - t0)./N;
I = length(n);

X = zeros(I, 1);	%final error
e = zeros(I, 1);	%global error

for i=[1:I]

	t=0;
	x = x0;
	dt = Dt(i);
	
	while t<tf - 0.5*dt;	
		x = RK3(x, t, dt, u);
		t = t+dt;
	end	
	
	X(i) = x;
	
end;

e = X - X(I);

plot(log(Dt), log(e), '-o');
Accuracy = polyfit(log(Dt(1:I-1)'), log(e(1:I-1)), 1)(1)

