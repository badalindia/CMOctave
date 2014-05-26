function [a0, a1] = Bound2(x0, x1, d0, d1)

% Given that x0 is bounded between 0 and d0 ------------ (1)
% Given that x1 is bounded between 0 and d1 ------------ (2)
% the function bounds (x1 - x0) between 0 and (d1 - d0) without breaking the conditions (1) and (2)

D = d1 - d0;

if (D >= 0)
	
	if( (x1 - x0) < 0 )
	
		%x0 = x1;
		m = (x0 + x1)/2;
		x0 = x1 = m;
	
	elseif ( (x1 - x0) > D )
	
		if ( d0 > 0 && d1 > 0)
		
			x1 = x0 + 0.99*D;
		
		elseif ( d0 < 0 && d1 < 0)
		
			x0 = x1 - 0.99*D;
		
		end;
		
		
	end;

elseif (D < 0)

	if( (x1 - x0) > 0 )
	
		%x1 = x0;
		m = (x0 + x1)/2;
		x0 = x1 = m;
	
	elseif ( (x1 - x0) < D )

		if ( d0 > 0 && d1 > 0)
		
			x0 = x1 - 0.99*D;
		
		elseif ( d0 < 0 && d1 < 0)
		
			x1 = x0 + 0.99*D;
		
		end;

	end;

end;

a0 = x0;
a1 = x1;
