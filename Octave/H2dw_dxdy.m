function [Hx, Hy] = H2dw_dxdy(x, y, f, fx, fy, fxy, x0, y0)

%w stands for warped

%note
% grid is assumed to be square  hx = hy = h
% grid is assumed to be periodic 

%Hermite basis f & g
Hf = @(x) (1 - 3.*x.^2 + 2.*x.^3);
Hg = @(x) x.*(1 - x).^2;
Hfx = @(x) (- 6.*x + 6.*x.^2);
Hgx = @(x) (1 + 3.*x.^2 - 4.*x);

%grid properties
[Ix, Iy] = size(x0);
h = x(1, 2) - x(1, 1);
L = length(x);
xDelta = x(1, L) - x(1, 1);
yDelta = y(L, 1) - y(1, 1);

INDEX = @(Ix, Iy) (Ix-1)*L + Iy;

I0x = floor((x0-x(1,1))./h) + 1;
I1x = I0x + 1;
I0y = floor((y0-y(1,1))./h) + 1;
I1y = I0y + 1;

	I0xPeriod = floor((I0x-1)/(L-1));
	I0xCell = mod(I0x-1, L-1) + 1;
	I1xPeriod = floor((I1x-1)/(L-1));
	I1xCell = mod(I1x-1, L-1) + 1;
	I0yPeriod = floor((I0y-1)/(L-1));
	I0yCell = mod(I0y-1, L-1) + 1;
	I1yPeriod = floor((I1y-1)/(L-1));
	I1yCell = mod(I1y-1, L-1) + 1;
	
%Ixy
I00 = INDEX(I0xCell, I0yCell);
I10 = INDEX(I1xCell, I0yCell);
I01 = INDEX(I0xCell, I1yCell);
I11 = INDEX(I1xCell, I1yCell);

dx = (x0 - (x(I00) + xDelta*I0xPeriod))./h;
dy = (y0 - (y(I00) + yDelta*I0yPeriod))./h;


%jump in 
fJumpX = f(:, L) - f(:, 1);
fJumpY = f(L, :)' - f(1, :)';

%fxy
f00 = f(I00)  + I0xPeriod.*fJumpX(I0yCell) + I0yPeriod.*fJumpY(I0xCell);
f01 = f(I01)  + I0xPeriod.*fJumpX(I1yCell) + I1yPeriod.*fJumpY(I0xCell);
f10 = f(I10)  + I1xPeriod.*fJumpX(I0yCell) + I0yPeriod.*fJumpY(I1xCell);
f11 = f(I11)  + I1xPeriod.*fJumpX(I1yCell) + I1yPeriod.*fJumpY(I1xCell);

i = [1:Ix*Iy];


Fx = [	Hfx(dx(i)) .* Hf(dy(i)); Hfx(dx(i)) .* Hf(1-dy(i)); Hfx(dx(i)) .* Hg(dy(i)); Hfx(dx(i)) .* -Hg(1-dy(i));
		-Hfx(1-dx(i)) .* Hf(dy(i)); -Hfx(1-dx(i)) .* Hf(1-dy(i)); -Hfx(1-dx(i)) .* Hg(dy(i)); -Hfx(1-dx(i)) .* -Hg(1-dy(i));
		Hgx(dx(i)) .* Hf(dy(i)); Hgx(dx(i)) .* Hf(1-dy(i)); Hgx(dx(i)) .* Hg(dy(i)); Hgx(dx(i)) .* -Hg(1-dy(i));
		Hgx(1-dx(i)) .* Hf(dy(i)); Hgx(1-dx(i)) .* Hf(1-dy(i)); Hgx(1-dx(i)) .* Hg(dy(i)); Hgx(1-dx(i)) .* -Hg(1-dy(i))]./h;
		
Fy = [	Hf(dx(i)) .* Hfx(dy(i)); 	Hf(dx(i)) .* -Hfx(1-dy(i)); 	Hf(dx(i)) .* Hgx(dy(i)); 	Hf(dx(i)) .* Hgx(1-dy(i));
		Hf(1-dx(i)) .* Hfx(dy(i)); 	Hf(1-dx(i)) .* -Hfx(1-dy(i)); 	Hf(1-dx(i)) .* Hgx(dy(i)); 	Hf(1-dx(i)) .* Hgx(1-dy(i));
		Hg(dx(i)) .* Hfx(dy(i)); 	Hg(dx(i)) .* -Hfx(1-dy(i)); 	Hg(dx(i)) .* Hgx(dy(i)); 	Hg(dx(i)) .* Hgx(1-dy(i));
		-Hg(1-dx(i)) .* Hfx(dy(i)); -Hg(1-dx(i)) .* -Hfx(1-dy(i)); 	-Hg(1-dx(i)) .* Hgx(dy(i)); -Hg(1-dx(i)) .* Hgx(1-dy(i))]./h;
		
		
		
M = [	      f00(i);	      f01(i);	  h*fy(I00)(i); 	  h*fy(I01)(i);
		      f10(i);	      f11(i);	  h*fy(I10)(i);		  h*fy(I11)(i);
		h*fx(I00)(i);	h*fx(I01)(i);	h*h*fxy(I00)(i);	h*h*fxy(I01)(i); 
		h*fx(I10)(i);	h*fx(I11)(i);	h*h*fxy(I10)(i);	h*h*fxy(I11)(i)];
		
		
Hx = zeros(Ix, Iy);
Hy = zeros(Ix, Iy);
Hx(i) = sum(Fx.*M)(i);
Hy(i) = sum(Fy.*M)(i);
