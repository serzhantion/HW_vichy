function [phi, dphi, ddphi] = coord_for_coloc(a,b,n,k)
syms 'x';
jac = myjacobi(k,n);
phi=sym(zeros(n,1));
dphi=sym(zeros(n,1));
ddphi=sym(zeros(n,1));
for i = 1:n
    phi(i) =(x-a)^2*(x-b)^2*jac(i);
    dphi(i) = diff(phi(i));
    ddphi(i) = diff(dphi(i));
end
end