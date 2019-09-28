function [Sx,Sy,fM] = fdm2d(f,p)
% model parameters
r = p.r;
rho = p.rho;
theta = p.theta;
T = p.T;
Nt = p.Nt;
sx = p.sigmax;
sy = p.sigmay;
ax = r - 0.5*sx^2;
ay = r-0.5*sy^2;
% discretization parameters
bx = p.bx;
Nx = p.Nx;
by = p.by;
Ny = p.Ny;
dx = 2*bx/(Nx-1);
xv = -bx:dx:bx;
dy = 2*by/(Nx-1);
yv = -by:dy:by;
dz = sqrt(dx^2+dy^2);
[f1,yD] = feval(f,xv,yv,p);
% derivatives over the boundaries
Dx0 = yD.Dx0;
DxN = yD.DxN;
Dy0 = yD.Dy0;
DyN = yD.DyN;
D00 = yD.D00;
D0N = yD.D0N;
DN0 = yD.DN0;
DNN = yD.DNN;
% combine columns
f0 = reshape(f1,Ny*Nx,1);
Sx = exp(xv);
Sy = exp(yv);
qp0 = 0.5*ax/dx + 0.5*sx^2/dx^2;
qm0 = -0.5*ax/dx + 0.5*sx^2/dx^2;
q0p = 0.5*ay/dx + 0.5*sy^2/dy^2;
q0m = -0.5*ax/dx + 0.5*sy^2/dy^2;
q00 = -r - sx^2/dx^2 - sy^2/dy^2;
qpp = 0.25*rho*sx*sy/dx/dy;
qmm=qpp;
qpm=-qpp;
qmp=-qpp;
I = ones(Nx,1);
I0 = ones(Nx-1,1);
% matrices for x-boundaries
D = spdiags(q00*I,0,Nx,Nx)+spdiags(qm0*I,-1,Nx,Nx)+spdiags(qp0*I,1,Nx,Nx);
C = spdiags(q0m*I,0,Nx,Nx)+spdiags(qmm*I,-1,Nx,Nx)+spdiags(qpm*I,1,Nx,Nx);
E = spdiags(q0p*I,0,Nx,Nx)+spdiags(qmp*I,-1,Nx,Nx)+spdiags(qpp*I,1,Nx,Nx);
% boundary corrections for x
D(1,2) = qp0+qm0;
D(Nx,Nx-1) = qp0+qm0;
C(1,2) = qpm+qmm;
C(Nx,Nx-1) = qpm+qmm;
E(1,2) = qpp+qmp;
E(Nx,Nx-1) = qpp+qmp;
% matricies for y-boundaries
B = E + spdiags(q0m*I,0,Nx,Nx)+spdiags(qmm*I,-1,Nx,Nx)+spdiags(qpm*I,1,Nx,Nx);
F = C + spdiags(q0p*I,0,Nx,Nx)+spdiags(qmp*I,-1,Nx,Nx)+spdiags(qpp*I,1,Nx,Nx);
% boundary corrections for y
B(1,2) = B(1,2) + qmm;
B(Nx,Nx-1) = B(Nx,Nx-1) + qpm;
F(1,2) = F(1,2) + qmp;
F(Nx,Nx-1) = F(Nx,Nx-1) + qpp;
% solve
Q = sparse(Nx*Ny,Nx*Ny);
G = sparse(Nx,Ny);
% first row
Q(1:Nx,1:Nx) = D;
Q(1:Nx,(Nx+1):(2*Nx)) = B;
G(:,1) = -2*(qmm+q0m+qpm)*Dy0*dy;
G(1,1) = G(1,1)-2*(qmm+qm0+qmp)*Dx0(1)*dx-2*qmm*D00*dz;
G(Nx,1) = G(Nx,1)+2*(qpm+qp0+qpp)*DxN(1)*dx-2*qpm*DN0*dx;
for jy=2:Ny-1
    Q((jy-1)*Nx+1:jy*Nx,(jy-1)*Nx+1:jy*Nx)=D;
    Q((jy-1)*Nx+1:jy*Nx,(jy-2)*Nx+1:(jy-1)*Nx)=C;
    Q((jy-1)*Nx+1:jy*Nx,jy*Nx+1:(jy+1)*Nx)=E;
    G(1,jy) = -2*(qmm+qm0+qmp)*Dx0(jy)*dx;
    G(Nx,jy) = 2*(qpm+qp0+qpp)*DxN(jy)*dx;
end
% last row
Q((Ny-1)*Nx+1:Ny*Nx,(Ny-1)*Nx+1:Ny*Nx) = D;
Q((Ny-1)*Nx+1:Ny*Nx,(Ny-2)*Nx+1:(Ny-1)*Nx) = F;
G(:,Ny) = 2*(qmm+q0m+qpm)*DyN*dy;
G(1,Ny) = G(1,Ny)-2*(qmm+qm0+qmp)*Dx0(Ny)*dx-2*qmp*D0N*dz;
G(Nx,Ny) = G(Nx,Ny)+2*(qpm+qp0+qpp)*DxN(Ny)*dx-2*qpp*DNN*dz;
g = reshape(G,Ny*Nx,1);
dt = T/(Nt-1);
QI = speye(Ny*Nx) - theta*dt*Q;
QE = speye(Ny*Nx) + (1-theta)*dt*Q;
f = f0;
for m=1:Nt-1
    f = QI\(QE*f + g*dt);
end
fM = reshape(f,Nx,Ny);




