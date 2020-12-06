function [t, y] = RK(t_0, h, t_n, y_0, A)%t_0 - первый узел, h - шаг, t_n - последний узел, y_0 - начальное условие, A - матрица
    t = t_0:h:t_n;
    y = zeros(rank(A), length(t));
    y(:,1) = y_0;
    for i = 2:length(t)
        k_1 = h*A*y(:,i-1);
        k_2 = h*A*(y(:,i-1)+k_1/2);
        k_3 = h*A*(y(:,i-1)+k_2/2);
        k_4 = h*A*(y(:,i-1)+k_3);
        y(:,i) = y(:,i-1) + (k_1 + 2*k_2 + 2*k_3 + k_4)/6;
    end
end
