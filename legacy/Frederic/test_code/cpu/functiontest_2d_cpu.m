xi = yi = -10.0
xf = yf = 10.0;
nx = 700;
ny = 700;
a = 2;
b = 2;
c = 1;
n = 1;
subplot(3,1,1);
mx = linspace(xi,xf,nx);
my = linspace(yi,yf,ny);
[xx,yy] = meshgrid(mx,my);
mz = cos(xx.^b)./(xx.^a+yy.^a+c).^n;
mesh(mx,my,mz)

subplot(3,1,2);
S=load("fort.2000");
plot3(S(:,1), S(:,2),S(:,3))
%mz = cos((xx+yy).^b) ./ (xx.^a+yy.^a+c).^n;
%mesh(mx,my,mz)

subplot(3,1,3);
quiver(S(:,1), S(:,2),S(:,4),S(:,5))