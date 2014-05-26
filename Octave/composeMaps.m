function [H, Hx, Hy, Hxy, V, Vx, Vy, Vxy] = composeMaps(X, Y, xDelta, yDelta, Xc, Yc, Hc, Hxc, Hyc, Hxyc, Vc, Vxc, Vyc, Vxyc, Xf, Yf, Hf, Hxf, Hyf, Hxyf, Vf, Vxf, Vyf, Vxyf)

eh = 10^-6;
[Iy Ix] = size(X);
i = j = [2:Ix-1];

%going through first diffeo
X0_xpyp = H2dw(Xc, Yc, Hc, Hxc, Hyc, Hxyc, X+eh, Y+eh);
X0_xpyn = H2dw(Xc, Yc, Hc, Hxc, Hyc, Hxyc, X+eh, Y-eh);
X0_xnyp = H2dw(Xc, Yc, Hc, Hxc, Hyc, Hxyc, X-eh, Y+eh);
X0_xnyn = H2dw(Xc, Yc, Hc, Hxc, Hyc, Hxyc, X-eh, Y-eh);

Y0_xpyp = H2dw(Xc, Yc, Vc, Vxc, Vyc, Vxyc, X+eh, Y+eh);
Y0_xpyn = H2dw(Xc, Yc, Vc, Vxc, Vyc, Vxyc, X+eh, Y-eh);
Y0_xnyp = H2dw(Xc, Yc, Vc, Vxc, Vyc, Vxyc, X-eh, Y+eh);
Y0_xnyn = H2dw(Xc, Yc, Vc, Vxc, Vyc, Vxyc, X-eh, Y-eh);

X1_xpyp = H2dw(Xf, Yf, Hf, Hxf, Hyf, Hxyf, X0_xpyp, Y0_xpyp);
X1_xpyn = H2dw(Xf, Yf, Hf, Hxf, Hyf, Hxyf, X0_xpyn, Y0_xpyn);
X1_xnyp = H2dw(Xf, Yf, Hf, Hxf, Hyf, Hxyf, X0_xnyp, Y0_xnyp);
X1_xnyn = H2dw(Xf, Yf, Hf, Hxf, Hyf, Hxyf, X0_xnyn, Y0_xnyn);

Y1_xpyp = H2dw(Xf, Yf, Vf, Vxf, Vyf, Vxyf, X0_xpyp, Y0_xpyp);
Y1_xpyn = H2dw(Xf, Yf, Vf, Vxf, Vyf, Vxyf, X0_xpyn, Y0_xpyn);
Y1_xnyp = H2dw(Xf, Yf, Vf, Vxf, Vyf, Vxyf, X0_xnyp, Y0_xnyp);
Y1_xnyn = H2dw(Xf, Yf, Vf, Vxf, Vyf, Vxyf, X0_xnyn, Y0_xnyn);

H = (X1_xpyp + X1_xnyn + X1_xpyn + X1_xnyp)/4;
Hx = (X1_xpyp - X1_xnyn + X1_xpyn - X1_xnyp)/(4*eh);
Hy = (X1_xpyp - X1_xnyn - X1_xpyn + X1_xnyp)/(4*eh);
Hxy = (X1_xpyp + X1_xnyn - X1_xpyn - X1_xnyp)/(4*eh*eh);

V = (Y1_xpyp + Y1_xnyn + Y1_xpyn + Y1_xnyp)/4;
Vx = (Y1_xpyp - Y1_xnyn + Y1_xpyn - Y1_xnyp)/(4*eh);
Vy = (Y1_xpyp - Y1_xnyn - Y1_xpyn + Y1_xnyp)/(4*eh);
Vxy = (Y1_xpyp + Y1_xnyn - Y1_xpyn - Y1_xnyp)/(4*eh*eh);
