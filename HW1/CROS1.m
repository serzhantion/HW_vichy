function [t, y] = CROS1(t_0, h, t_n, y_0, A)%t_0 - первый узел, h - шаг, t_n - последний узел, y_0 - начальное условие, A - матрица
    s = rank(A);
    E = eye(size(A));
    t = t_0:h:t_n;
    y = zeros(s, length(t));
    y(:,1) = y_0;
    for i = 2:length(t)
        w = (E-h*A*(1+1i)/2)\(A*y(:,i-1));%решение СЛАУ
        y(:,i) = y(:,i-1) + h*(w+conj(w))/2;%получение конечного ответа
    end
end
