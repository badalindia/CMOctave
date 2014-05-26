function [H, Hx] = H1d(x, f, fx, x0)

% 1D Hermite Interpolation
% x is the grid
% f is the function values at those grid points
% fx  = df/dx
% the function returns the interpolated value at x0
	%x0 should lie within the domain of x


xN = size(x)(2);

%Hermite interpolant basis
Hf = @(x) (1 - 3.*x.^2 + 2.*x.^3);
Hg = @(x) x.*(1 - x).^2;
Hfx = @(x) - 6.*x + 6.*x.^2;
Hgx = @(x) (1 - x).^2 - 2.*x.*(1-x);

h_ = x(2) - x(1);

%problem parameters
I0 = floor((x0-x(1))./h_) + 1;
I1 = I0 + 1;

dx = (x0 - x(I0))./h_;

%Hermite interpolant
H = f(I0).*Hf(dx) + f(I1).*Hf(1 - dx) + h_*(fx(I0).*Hg(dx) - fx(I1).*Hg(1 - dx));
Hx = (f(I0).*Hfx(dx) - f(I1).*Hfx(1 - dx))/h_ + fx(I0).*Hgx(dx) + fx(I1).*Hgx(1 - dx);


end


