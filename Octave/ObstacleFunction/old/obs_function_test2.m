N = 128;
x = linspace(0,1,N+1);
[X, Y] = meshgrid(x);
h = x(2) - x(1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#buffering matrix for implicite source term calculations
[Pc00, Pc01, Pc10, Pc11] = obs_funct_implicite2(X, Y, h);
