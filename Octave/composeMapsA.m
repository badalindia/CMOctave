function [H, Hx, Hy, Hxy, V, Vx, Vy, Vxy] = composeMapsA(X, Y, Xc, Yc, Hc, Hxc, Hyc, Hxyc, Vc, Vxc, Vyc, Vxyc, Xf, Yf, Hf, Hxf, Hyf, Hxyf, Vf, Vxf, Vyf, Vxyf)

%A stands for 'Analytically'


H = X.*0;
Hx = X.*0;
Hy = X.*0;
Hxy = X.*0;
V = X.*0;
Vx = X.*0;
Vy = X.*0;
Vxy = X.*0;


%passing through map 1
[h1, h1x, h1y, h1xy] = H2dwF(Xc, Yc, Hc, Hxc, Hyc, Hxyc, X, Y);
[v1, v1x, v1y, v1xy] = H2dwF(Xc, Yc, Vc, Vxc, Vyc, Vxyc, X, Y);

[h2, h2x, h2y, h2xy, h2xx, h2yy] = H2dwF2(Xf, Yf, Hf, Hxf, Hyf, Hxyf, h1, v1);
[v2, v2x, v2y, v2xy, v2xx, v2yy] = H2dwF2(Xf, Yf, Vf, Vxf, Vyf, Vxyf, h1, v1);

H = h2;
Hx = h2x.*h1x + h2y.*v1x;
Hy = h2x.*h1y + h2y.*v1y;
Hxy = h2xx.*h1x.*h1y + h2xy.*(h1x.*v1y + h1y.*v1x) + h2yy.*v1x.*v1y + h2x.*h1xy + h2y.*v1xy;

V = v2;
Vx = v2x.*h1x + v2y.*v1x;
Vy = v2x.*h1y + v2y.*v1y;
Vxy = v2xx.*h1x.*h1y + v2xy.*(h1x.*v1y + h1y.*v1x) + v2yy.*v1x.*v1y + v2x.*h1xy + v2y.*v1xy;
