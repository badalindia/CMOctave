function _G = bilinearInterpolation(x, y, _x, _y, G)

[Iy Ix] = size(x);
xMin = x(1, 1);
xMax = x(1, Ix);
xDelta = xMax - xMin;
yMin = y(1, 1);
yMax = y(Iy, 1);
yDelta = yMax - yMin;
h = x(1, 2) - x(1, 1);

		L = length(x);

		I0x = floor((_x - xMin)./h) + 1;
		I1x = I0x + 1;
		I0y = floor((_y - yMin)./h) + 1;
		I1y = I0y + 1;
		
		ifI0xW = (I0x<1) - (I0x>L);
		ifI1xW = (I1x<1) - (I1x>L);
		ifI0yW = (I0y<1) - (I0y>L);
		ifI1yW = (I1y<1) - (I1y>L);
		
		I0x = I0x + (L-1)*ifI0xW;
		I1x = I1x + (L-1)*ifI1xW;
		I0y = I0y + (L-1)*ifI0yW;
		I1y = I1y + (L-1)*ifI1yW;
		
		I00 = (I0x-1)*L + I0y;	%Ixy
		I01 = (I1x-1)*L + I0y;
		I10 = (I0x-1)*L + I1y;
		I11 = (I1x-1)*L + I1y;
		
		dx = (_x - (x(I00) - (ifI0xW)*(x(1,L) - x(1,1))))./h;
		dy = (_y - (y(I00) - (ifI0yW)*(y(L,1) - y(1,1))))./h;
	
		_Ixy = [1:Ix*Iy];
		_G = (G(I00).*(1-dx) + G(I01).*(dx)).*(1-dy) + (G(I10).*(1-dx) + G(I11).*(dx)).*(dy);
		
