

Nc = 32;
Nf = 64;
Np = 64;

x = linspace(0, 1, Np+1);

[X, Y] = meshgrid(x, x);

C = 2;
Q = @(X, Y) sin(C*pi*X).*sin(C*pi*Y);				

h = x(2) - x(1);

U = @(X, Y, T) -cos(T).*(sin(pi.*X).^2).*sin(2.*pi.*Y);
V = @(X, Y, T)  cos(T).*(sin(pi.*Y).^2).*sin(2.*pi.*X);

Diffeo2D(X, Y, Q, U, V, Nc, Nf, 0, 2*pi, h/2, true);



