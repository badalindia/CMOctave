function ans = Bound(x, b)

%bound x between 0 and b

if(b >= 0)

	if(x < 0)
		x = 0;
	end;
	
	if(x > b)
		x = 0.999*b;
	end;

end;

if(b < 0)

	if(x > 0)
		x = 0;
	end;
	
	if(x < b)
		x = 0.999*b;
	end;

end;

ans = x;
