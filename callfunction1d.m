function [y,by,dy] = callfunction1d(x,p)
K = p.K;
dx = x(2)-x(1);
y = max(exp(x) - K,0);
by = NaN; %Dirichlet
dy = [0,exp(x(end))]; %Neumann