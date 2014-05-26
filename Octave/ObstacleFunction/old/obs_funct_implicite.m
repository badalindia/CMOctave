%onstacle function 
function [Pc00, Pc01, Pc10, Pc11] = obs_funct_implicite(x, y, dx)

%[x, y] is the grid

x0=0.5;
y0=0.5;
r0=0.15;

phi	 = -(sqrt((x-x0).^2+(y-y0).^2)-r0);
phix = (x-x0)./sqrt((x-x0).^2+(y-y0).^2 + 0.001);
phiy = (y-y0)./sqrt((x-x0).^2+(y-y0).^2 + 0.001);

eps=2*dx;
Ch=(tanh(phi./eps)+1).*0.5;
M = (tanh((r0/2-phi)./(eps*4))+1).*0.5;

#{
contourf(x, y, M, 50)
axis([0, 1, 0, 1], "square");
shading flat
colorbar;
drawnow;
print('M.png', '-r200');
#}

Pc00 = Ch.*(1-M.*phiy.*phiy);
Pc01 = Pc10 = Ch.*(M.*phix.*phiy);
Pc11 = Ch.*(1-M.*phix.*phix);




