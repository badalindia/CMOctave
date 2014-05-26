function [xn, yn] = AdvectionUV_Euler1D(Xc, Yc, x, y, dt, U, Ux, Uy, Uxy, V, Vx, Vy, Vxy)

Un = H2dw(Xc, Yc, U, Ux, Uy, Uxy, x, y);
Vn = H2dw(Xc, Yc, V, Vx, Vy, Vxy, x, y);

xn = x - dt*Un;
yn = y - dt*Vn;
