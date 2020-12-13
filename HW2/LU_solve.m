function x = LU_solve(A, b)%A - матрица, b - вектор
[L, U] = LU(A);
s = size(A);
r = s(1);
x = zeros(1,r);
y = zeros(1,r);
y(1) = b(1);
%решение Ly = b
for i = 2:r
    y(i) = b(i);
    for j =1:i-1
        y(i) = y(i) - L(i,j)*y(j);
    end
end
%решение Ux = y
x(r) = y(r)/U(r,r);
for i = r-1:-1:1
    x(i) = y(i);
    for j = r:-1:i+1
        x(i) = x(i) - U(i,j)*x(j);
    end
    x(i) = x(i)/U(i,i);
end
x = transpose(x);
end
