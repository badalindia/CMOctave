Np = 128;
x = linspace(0, 1, Np);
[X, Y] = meshgrid(x);

f = exp(-5000*(((X-0.5).^2 + (Y-0.5).^2) - 0.15).^2);
f = f.*(X < 0.5);


max(max(f));
min(min(f));

%surf(X, Y, f);


contour(X, Y, f, [0.90, 0.90]);
axis([0, 1, 0, 1], 'square');

hold on;
plot(X, X', 'r', "linewidth", 1);
plot(Y, Y', 'r', "linewidth", 1);
hold off;
		
print('test.png', '-r400');


