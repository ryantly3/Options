p.r = 0.04;
p.rho = 0.30;
p.sigmax = 0.20;
p.sigmay = 0.30;
p.K = 1.00;
p.theta = 0.50;
p.T = 0.50;
p.Nt = 30;
p.bx = 0.50;
p.Nx = 51;
p.by = 0.50;
p.Ny = 51;
[Sx,Sy,fM] = fdm2d(@callfunction2d,p);
ix = find((Sx>0.8).*(Sy<1.2));
iy = find((Sy>0.8).*(Sy<1.2));
contour(Sx(ix),Sy(iy),fM(iy,ix)');