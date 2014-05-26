function temp = T(x, y, coeff)

dx=x(1,2)-x(1,1);
dy=dx;

x0=0.5;
y0=0.125;

sigma1=1/15;


sigma2x=sigma1*2.5*coeff;
sigma2y=sigma1*2.5/coeff;

rho=x.*0;
temp=rho;

for i=-1:1
    for j=-1:1
        rho= rho + exp(-(((x-x0+i)/sigma1).^8+((y-y0+j)/sigma1).^8)) + exp(-(((x-x0+i)/sigma1/2).^4+((y-y0+j)/sigma1/2).^4));
        temp = temp + exp(-(((x-x0+i)/sigma1).^8+((y-y0+j)/sigma1).^8)) + exp(-(((x-x0+i)/(sigma2x)).^2+((y-y0+j)/sigma2y).^2));
    end
end




