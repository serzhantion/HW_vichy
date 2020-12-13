function x = QR_solve(A,b)%A - матрица, b - вектор
[Q, R] = QR(A);
s = size(A);
r = s(1);
b = transpose(Q)*b;
%решение Rx = Q^T * b
x(r) = b(r)/R(r,r);
for i = r-1:-1:1
    x(i) = b(i);
    for j = r:-1:i+1
        x(i) = x(i) - R(i,j)*x(j);
    end
    x(i) = x(i)/R(i,i);
end
x = transpose(x);
end
