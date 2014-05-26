function myquiver(n, x, y, u, v,path)

figure(n);
quiver(x, y, u, v);
axis([min(min(x)), max(max(x)), min(min(y)), max(max(y))], "square");
drawnow;
print(path, '-r200');


