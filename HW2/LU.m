function [L, U] = LU(A)%А - матрица
s = size(A);
n = s(1);
L = eye(s);
U = A;
for i = 1:n-1
    M = diag(diag(eye(s)));
    M(i+1:n,i) = -A(i+1:n,i)/A(i,i);
    L = L*M^(-1);
    U = M*U;
end
end
