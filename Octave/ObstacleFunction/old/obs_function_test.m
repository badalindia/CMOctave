N = 128;

[x, y] = meshgrid(linspace(0, 1, N+1));
h = 1/N;

dt = 10^-2;
eta = 10^-3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#buffering matrix for implicite source term calculations
[Pc00, Pc01, Pc10, Pc11] = obs_funct_implicite(x, y, h);

%adding 1/dt and inverting
P00_ = 1/dt+1/eta*Pc00;
P01_ =      1/eta*Pc01;
P10_ =      1/eta*Pc10;
P11_ = 1/dt+1/eta*Pc11;

D_ = P00_.*P11_ - P01_.*P01_;

P00 =  P11_./D_;
P01 = -P01_./D_;
P10 = -P10_./D_;
P11 =  P00_./D_;

u = v = x.*0+1;

u2 = 1/dt*(P00.*u + P01.*v);
v2 = 1/dt*(P10.*u + P11.*v);

th = linspace(0,2*pi,100)'; 
circsx = 0.15.*cos(th) + 0.5; 
circsy = 0.15.*sin(th) + 0.5; 

quiver(x, y, u2, v2);
axis([0, 1, 0, 1], "square");
hold on;
plot(circsx, circsy);
hold off;
print('v(eta=10).png', '-r400');
