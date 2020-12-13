function [x] = zeidel(A, b, epsilon)
%A - матрица, b - вектор, epsilon - заданная погрешность
    disp('Метод Зейделя');
    s = size(A);
    E = eye(s);
    D = diag(diag(A));
    H = eye(s) - inv(D) * A;
    g = inv(D)*b;
    H_l = tril(H, -1);
    H_r = triu(H);
    % приведение к виду метода простой итерации
    H = inv(E - H_l)*H_r;
    g = inv(E - H_l)*g;
    x_0 = zeros(s(1), 1);
    i = 0;
    if (max(max(abs(H))) <= 1)%достаточное условие сходимости
        x = H*x_0 + g;
        i = i +1;
        while norm(x - x_0) > epsilon
            x_0 = x;
            x = H * x + g;
            i = i + 1;
        end
    else
          disp('Ошибка: спектральный радиус > 1');
          x = x_0;
    end
    disp('Количество итераций: ')
    disp(i);
end
