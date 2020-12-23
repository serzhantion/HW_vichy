function y = ritz(a,b,alpha, beta, n, p, q, f)
k = 1;
%test2
syms 'fun(x)';
fun(x) = x;
phi=sym(zeros(n,1));
dphi=sym(zeros(n,1));
phi(1) = 1;
phi(2) = fun;
dphi(1) = 0;
dphi(2) = diff(fun);
[phi(3:n), dphi(3:n), ~]=myphi_dphi_d2kk_good_anal(k,n-2);

%test1
%[phi, dphi, ~]=myphi_dphi_d2kk_good_anal(k,n);
d = zeros(n,1);
A = zeros(n);

for i = 1:n
    for j = 1:n
        [q_r, q_l] = bc(alpha,beta, a, b, p, phi(i), phi(j));
        A(i,j) = double(int(p*dphi(i)*dphi(j) + q*phi(i)*phi(j), a, b)) + q_r + q_l;
    end
    func = f*phi(i);
    d(i) = int(func, [a b]);
end
digits(4);
alpha = vpa(A\d);
y = alpha.*phi;
y = sum(y);
end