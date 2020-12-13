function [Q,R] = QR(A)%A - матрица
s = size(A);
n = s(1);
Q = eye(s);
R = A;
for j = 1:n-1
    for i = n:-1:j+1
        %задание матрицы поворотов
        T = eye(s);
        theta = atan(-R(i,j)/R(j,j));
        T(i,i) = cos(theta);
        T(j,j) = T(i,i);
        T(i,j) =  sin(theta);
        T(j,i) =  -T(i,j);
        R = T*R;% => R(i,j)=0
        Q = Q*transpose(T);
    end
end
end
