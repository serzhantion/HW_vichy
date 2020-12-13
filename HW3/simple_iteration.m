function [x] = simple_iteration(A, b, epsilon)
%A - матрица, b - вектор, epsilon - заданная погрешность
    disp('Метод простой итерации')
    s = size(A);
    D = diag(diag(A));
    % приведение к виду x_k-1 = H*x_k+g
    H = eye(s) - inv(D) * A;
    g = inv(D)*b;
    x_0 = zeros(s(1), 1);
    x = ones(s(1),1);
    i = 0;
    if (max(abs(eig(H))) < 1)%н. и д. условие сходимости
        x = H*x_0 + g;
        i = i +1;
        while norm(x - x_0) > epsilon
            x_0  = x;
            x = H * x + g;
            i = i + 1;
        end
    else
        disp('Ошибка: спектральный радиус > 1');
        x = x_0;
    end
    disp('Количество итераций:');
    disp(i);
end
