clear;
clc;
p.theta = 0.50;
p.r = 0.04;
p.sigma = 0.30;
p.t = 0.25;
p.K = 0.95;
p.tnumber = 60;
p.xboundary = 0.40;
p.xnumber = 80;
p.boundtype = 0;
[xv,tv,FT,exact,error] = fdm1d(@putfunction1d,p,2);
surf(exp(xv),tv,FT');