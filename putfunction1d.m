function [y,by,dy] = putfunction1d(x,p)
K = p.K;
dx = x(2)-x(1);
y = max(K-exp(x),0);
by = NaN; %Dirichlet
dy = [-exp(x(1)),0]; %Neumann