clear all;
nrows = 128;


S2d = load('rho_16_9.blk', '-ascii');

size(S2d)

S2d = reshape(S2d, nrows+1, nrows+1);

i = [1:nrows];
j = [1:nrows];

for k = [1:nrows]
	S(k, i, j) = S2d(i, j);
end;


[min(min(min(S))) max(max(max(S)))]
data=permute(reshape(S,nrows,nrows,nrows),[1 2 3]);



X=linspace(0,1,nrows);
[x,y,z] = meshgrid(X,X,X);

figure;
p = patch(isosurface(x,y,z,data,1));
isonormals(x,y,z,data,p)
set(p,'FaceColor','red','EdgeColor','none');
daspect([1 1 1])
axis([0 1 0 1 0 1])
view(3);
camlight 
lighting gouraud
box on

