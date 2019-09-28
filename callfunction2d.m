function [y,yD] = callfunction2d(xv,yv,p)
K = p.K;
[sx,sy] = meshgrid(exp(xv),exp(yv));
dx = xv(2)-xv(1);
dy = yv(2)-yv(1);
y = max(min(sx,sy)-K,0);
%partial derivatives
yD.Dx0 = (y(2,:)-y(1,:))/dy;
yD.DxN = (y(end,:)-y(end-1,:))/dy;
yD.Dy0 = (y(:,2)-y(:,1))/dx;
yD.DyN = (y(:,end)-y(:,end-1))/dx;
%diagonal derivatives
dz = sqrt(dx^2+dy^2);
yD.D00 = (y(1,1)-y(2,2))/dz;
yD.D0N = (y(1,end)-y(2,end-1))/dz;
yD.DN0 = (y(end,1)-y(end-1,2))/dz;
yD.DNN = (y(end,end)-y(end-1,end-1))/dz;