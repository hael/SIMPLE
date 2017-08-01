xi = yi = -10.0
xf = yf = 10.0;
subplot(3,1,1);
nx = 1000;
ny = 1000;
a = 2;
b = 2;
c = 1;
n = 1;
mx = linspace(xi,xf,nx);
my = linspace(yi,yf,ny);
[xx,yy] = meshgrid(mx,my);
mz = cos(xx.^b)./(xx.^a+yy.^a+c).^n + cos((xx-5).^b)./((xx-5).^a+(yy-5).^a+c).^n+ cos((xx+5).^b)./((xx+5).^a+(yy+5).^a+c).^n;
mesh(mx,my,mz)

subplot(3,1,2);
S=load("fort.2000");
plot3(S(:,1), S(:,2),S(:,3))
%mz = cos((xx+yy).^b) ./ (xx.^a+yy.^a+c).^n;
%mesh(mx,my,mz)

subplot(3,1,3);
%colormap('default');
%mz = cos((xx+yy).^b) ./ (xx.^a+yy.^a+c).^n;
mz = cos((xx+yy).^b) ./ (xx.^a+yy.^a+c).^n;
mesh(mx,my,mz)

%subplot(3,1,3);
%S=load("fort.3000");
%plot3(S(:,1), S(:,2),S(:,3))
