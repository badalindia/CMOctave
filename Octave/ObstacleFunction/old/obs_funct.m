%onstacle function 
function [Fx, Fy] = obs_funct(X, Y, x, y, f, fx, fy, fxy)

%[X, Y] is the grid
%(x, y) are the points at which obstacle function value is calculated

dx=X(1, 2)-X(1, 1);

x0=0.5;
y0=0.5;
r0=0.15;

phi	= -(sqrt((x-x0).^2+(y-y0).^2)-r0);
phix = (x-x0)./sqrt((x-x0).^2+(y-y0).^2);
phiy = (y-y0)./sqrt((x-x0).^2+(y-y0).^2);

eps=dx*2;
Ch=(tanh(phi./eps)+1).*0.5;
M = (tanh((r0/2-phi)./(eps*8))+1).*0.5;

eh = 10^-6;
fx2 = (H2d_Warped(X, Y, f, fx, fy, fxy, x+eh, y) - H2d_Warped(X, Y, f, fx, fy, fxy, x-eh, y))./(2*eh);
fy2 = (H2d_Warped(X, Y, f, fx, fy, fxy, x, y+eh) - H2d_Warped(X, Y, f, fx, fy, fxy, x, y-eh))./(2*eh);
u = -fy2;
v =  fx2;

%un = (u.*phix.*phix + v.*phix.*phiy);
%vn = (u.*phix.*phiy + v.*phiy.*phiy);

ut =   u.*phiy.*phiy - v.*phix.*phiy;
vt = - u.*phix.*phiy + v.*phix.*phix;

un = u - ut.*M;
vn = v - vt.*M;

Fx = -Ch.*un;
Fy = -Ch.*vn;



