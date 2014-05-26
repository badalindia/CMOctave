%penalty function for Navier Stokes
function P = obs_functNS(x, y, dx)	

%[x, y] is the grid

x0=0.5;
y0=0.5;
r0=0.15 - dx;

phi	 = -(sqrt((x-x0).^2+(y-y0).^2)-r0);

eps=1.5*dx;
P = (phi>0); 
P = (abs(phi)<=eps).*0.5.*(1+phi./eps+(1/pi)*sin(pi.*phi./eps)) + (abs(phi)>eps).*P;





