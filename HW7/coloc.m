function y = coloc(a,b,n, p, q, f,v)
ch = zeros(1,n);
for i=1:n
        ch(i)=(a+b)/2 + (a-b)/2*cos((2*i-1)*pi/(n*2));
end
k = 1;
%test2
% syms 'x';
% phi=sym(zeros(n,1));
% dphi=sym(zeros(n,1));
% ddphi=sym(zeros(n,1));
% phi(1) = 7/2*x^3 - 17/2*x + x^2;
% phi(2) = -1/7*x^2 -2/7*x + 1;
% dphi(1) = diff(phi(1));
% dphi(2) = diff(phi(2));
% ddphi(1) = diff(dphi(1));
% ddphi(2) = diff(dphi(2));
% [phi(3:n), dphi(3:n), ddphi(3:n)]=coord_for_coloc(a,b,n-2,k);

%test1
[phi, dphi, ddphi]=myphi_dphi_d2kk_good_anal(k,n);
d = zeros(n,1);
A = zeros(n);
for i = 1:n
    for j = 1:n
        func = -p*ddphi(j) - diff(p)*dphi(j)+v*dphi(j) + q*phi(j);
        A(i,j) = func(ch(i));
    end
    d(i) = f(ch(i));
end
digits(4);
alpha = vpa(A\d);
y = alpha.*phi;
y = sum(y);
end