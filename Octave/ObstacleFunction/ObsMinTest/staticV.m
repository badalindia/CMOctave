N = 64;
x = linspace(0, 1, N+1);
h = x(2) - x(1);
[X, Y] = meshgrid(x);

%initialization of velocity field
u = sin(pi*X).^2;
v = sin(pi*Y).^2;

%obstacle function
[Pc00, Pc01, Pc10, Pc11] = obs_funct2(X, Y, h);

dt = 10^-3;
eta = 10^-2;

t=0;
tf = 1;
while t<tf

	uNew = u - dt/eta*(Pc00.*u + Pc01.*v);
	vNew = v - dt/eta*(Pc10.*u + Pc11.*v);
	
	u = uNew;
	v = vNew;
	
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
