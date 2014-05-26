N = 128;

[x, y] = meshgrid(linspace(0, 1, N+1));
h = 1/N;

dt = 10^-2;
eta = 10^-3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#buffering matrix for implicite source term calculations
[Pc00, Pc01, Pc10, Pc11] = obs_funct2(x, y, h);

u = v = x.*0+1;

u2 = dt/eta*(Pc00.*u + Pc01.*v);
v2 = dt/eta*(Pc10.*u + Pc11.*v);

th = linspace(0,2*pi,100)'; 
circsx = 0.15.*cos(th) + 0.5; 
circsy = 0.15.*sin(th) + 0.5; 

quiver(x, y, u2, v2);
axis([0, 1, 0, 1], "square");
hold on;
plot(circsx, circsy);
hold off;
print('v(eta=10).png', '-r400');
