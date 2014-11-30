function [H, Hx] = refineMap(x, f, fx, flag);

	%%%%%%%%% caution!!! this function assumes that the end point of the array is same as the first point

%overshooting possibility
L = size(x)(2);
h = x(2) - x(1);

i = [1+flag:2:L-1];

	%alternate neihbours

in_jump = (i-2)<1;
in = i-2 + in_jump*(L-1);
ip_jump = (i+2)>(L-1);
ip = i+2 - ip_jump*(L-1);

J = h*(L-1);

fi = f(i);
fin = f(in) - in_jump*J;
fip = f(ip) + ip_jump*J;


		%next neighbours
		i2n_jump = (i-1)<1;
		i2n = i-1 + i2n_jump*(L-1);
		
		i2p_jump = (i+1)>(L-1);
		i2p = i+1 - i2p_jump*(L-1);
		
		fi2n = f(i2n) - i2n_jump*J;
		fi2p = f(i2p) + i2p_jump*J;
		fxi2n = fx(i2n);
		fxi2p = fx(i2p);


df1 = abs(fip - fi)./(2*h);
df2 = abs(fi - fin)./(2*h);
df3 = abs(fi2p - fi)./h;
df4 = abs(fi - fi2n)./h;

df = max([df1, df2]);

shoot = (abs(fx(i)) > 3*df)*1;
%shoot2 = (fx(i)./(fip - fi) < 0) + (fx(i)./(fi - fin) < 0);
shoot2 = (fx(i)<0).*1;

if(sum(shoot2) > 0)
	i
	shoot2
	su
end;

	
	


H = f; Hx = fx;

H(L) = H(1)+h*(L-1);
Hx(L) = Hx(1);

if(sum(shoot)>0)

	
	for j=[1:length(i)]
		if(shoot(j))
			[H(i(j)), Hx(i(j))] = H1dw( [0, 2*h] , [f(i2n(j)), f(i2p(j))], [fx(i2n(j)), fx(i2p(j))], h);
		end;
	end
end

