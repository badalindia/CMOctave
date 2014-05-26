%onstacle function 
function [Pc00, Pc01, Pc10, Pc11] = obs_funct_implicite2(x, y, dx)

%[x, y] is the grid

x0=0.5;
y0=0.5;
r0=0.15;

phi	 = -(sqrt((x-x0).^2+(y-y0).^2)-r0);
phix = (x-x0)./sqrt((x-x0).^2+(y-y0).^2 + 0.001);
phiy = (y-y0)./sqrt((x-x0).^2+(y-y0).^2 + 0.001);

eps=3*dx;
Ch = (phi>0); 
Ch = (abs(phi)<=eps).*0.5.*(1+phi./eps+(1/pi)*sin(pi.*phi./eps)) + (abs(phi)>eps).*Ch;

phi2 = r0/2-phi;
eps2=4*dx;
M = (phi2>0); 
M = (abs(phi2)<=eps2).*0.5.*(1+phi2./eps2+(1/pi)*sin(pi.*phi2./eps2)) + (abs(phi2)>eps2).*M;

#{
contourf(x, y, M, 50)
axis([0, 1, 0, 1], "square");
shading flat;
colorbar;
drawnow;
print('M.png', '-r200');
#}

Pc00 = Ch.*(1-M.*phiy.*phiy);
Pc01 = Pc10 = Ch.*(M.*phix.*phiy);
Pc11 = Ch.*(1-M.*phix.*phix);




