function [xv, tv, FT, exact, error] = fdm1d(f,p,option)
theta = p.theta;
r = p.r; % risk free rate
sigma = p.sigma; % volatility
a = r - 0.5*sigma*sigma;
T = p.t; % maturity
Nt = p.tnumber; % time interval
dt = T/Nt;
tv = [0:dt:T];
bx = p.xboundary;
Nx = p.xnumber; % log-price interval
dx = 2*bx/Nx;
xv = [-bx:dx:bx];
boundtype = p.boundtype; %boundary condition
[f0,bf0,df0] = feval(f,xv',p);
qp = 0.5*a/dx + 0.5*sigma*sigma/dx/dx;
qm = -0.5*a/dx + 0.5*sigma*sigma/dx/dx;
q0 = -r - sigma*sigma/dx/dx;
%tridiagonal matrix Q
Q = diag(q0*ones(Nx+1,1)) + diag(qp*ones(Nx,1),1) + diag(qm*ones(Nx,1),-1);
g = zeros(Nx+1,1);
%boundary conditions
if boundtype % Dirichlet conditions
    g(1) = qm*bf0(1);
    g(end) = qp*bf0(2);
else % Neumann conditions
    Q(1,2) = qm + qp;
    Q(Nx+1,Nx) = qm + qp;
    g(1) = -2*dx*qm*df0(1);
    g(end) = 2*dx*qp*df0(2);
end
FT = zeros(Nx+1,Nt+1); % results matrix
call = zeros(Nx+1,Nt+1);
put = zeros(Nx+1,Nt+1);
FT(:,1) = f0; % initial condition
[call(:,1), put(:,1)] = blsprice(exp(xv),0.95,r,0,sigma);
for tndx = 2:Nt+1
    FT(:,tndx) = (eye(Nx+1)-theta*Q*dt)\((eye(Nx+1) + (1-theta)*Q*dt)*FT(:,tndx-1)+g*dt);
    [call(:,tndx), put(:,tndx)] = blsprice(exp(xv),0.95,r,tndx*dt,sigma); % exact solution
end
% error
call_error = norm(call - FT)/norm(call);
put_error = norm(put - FT)/norm(put);
% output
if option == 1
    exact = call;
    error = call_error;
else
    exact = put;
    error = put_error;
end