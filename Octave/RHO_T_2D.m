
x = linspace(0, 1, 129);
dy = dx = x(2) - x(1);
[x, y] = meshgrid(x);

x0=0.5;
y0=0.125;

sigma1=1/15;
coeff = 1.5;

sigma2x=sigma1*2.5*coeff;
sigma2y=sigma1*2.5/coeff;

rho=x.*0;
temp=rho;

for i=-1:1
    for j=-1:1
        rho= rho + exp(-(((x-x0+i)/sigma1).^8+((y-y0+j)/sigma1).^8)) + exp(-(((x-x0+i)/sigma1/2).^4+((y-y0+j)/sigma1/2).^4));
        temp=temp+ (exp(-(((x-x0+i)/sigma1).^8+((y-y0+j)/sigma1).^8)) + exp(-(((x-x0+i)/(sigma2x)).^2+((y-y0+j)/sigma2y).^2)));
    end
end


contourf(x, y, -rho, 50);
axis([0, 1, 0, 1], 'square');
shading flat;
colormap(gray);




