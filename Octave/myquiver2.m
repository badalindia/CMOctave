function myquiver2(n, x, y, u, v,path)

x0=0.55;
y0=0.5;
r0=0.15;

th = linspace(0,2*pi,100)'; 
circsx = r0.*cos(th) + x0; 
circsy = r0.*sin(th) + y0; 

dx = 1/64;
eps=1.5*dx;
circsx2 = (r0+eps).*cos(th) + x0; 
circsy2 = (r0+eps).*sin(th) + y0; 

circsx3 = (r0-eps).*cos(th) + x0; 
circsy3 = (r0-eps).*sin(th) + y0; 

figure(n);
quiver(x, y, u, v, 2);
axis([min(min(x)), max(max(x)), min(min(y)), max(max(y))], "square");
	hold on;
	plot(circsx, circsy, 'r', 'linewidth', 2);
	plot(circsx2, circsy2, 'b', 'linewidth', 1);
	plot(circsx3, circsy3, 'b', 'linewidth', 1);
	drawnow;
	hold off;
print(path, '-r400');


