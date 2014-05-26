function myplot(n, x, y, f, nC, map, range, path, res)

figure(n);
contourf(x, y, f, nC);
%pcolor(x, y, f);
colormap(map);

if(length(range)>1)
	caxis(range)
end

axis([min(min(x)), max(max(x)), min(min(y)), max(max(y))], "square");
shading flat;
colorbar;
drawnow;
print(path, res);


