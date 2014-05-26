function [xn, yn] = Advection_Euler1D(Xc, Yc, x, y, dt, U, V)

xn = x - dt*U;
yn = y - dt*V;
