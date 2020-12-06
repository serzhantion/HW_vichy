function [y, t, d] = Richardson(Number,met, t_0, t_n, y_0, A)

%Number - номер файла output(Number).dat, met - название используемого метода
%t_0 - первый узел, t_n - последний узел, y_0 - начальное условие, A - матрица
    
    if met == "Adams"%вычисление начального шага, в зависимости от метода
        h = 1e-4;
        por_met = 3;
    elseif met == "RK"
        h = 2.78/(max(abs(eig(A))));
        por_met = 4;
    elseif met == "CROS1"
        h = 1e-3;
        por_met = 2;
    end
    epsilon = 1e-10;%точность, до которой будем вычислять
    %применение метода для получения значения на первой сетке
    [t_1, v_1] = feval(met,t_0, h, t_n, y_0, A);
    d = ones(1,rank(A));
    file = fopen(strcat('output', int2str(Number), '.dat'),'w+');%откритие файла для записи
    
    %условие точности и дополнительно условие на количество узлов
    while and(length(t_1)<1e7, max(abs(d))>epsilon)
        h = h/2;%сгущаем сетку
        [t_2, v_2] = feval(met,t_0, h, t_n, y_0, A);%вычисляем значения на новой сетке
        dd = zeros(3,rank(A));%матрица с тремя векторами погрешностей
        
        %нахождение максимальной погрешности между сетками
        dd(1,:) = (v_2(:,1) - v_1(:,1))/(2^por_met - 1);
        m = 0;
        for j = 2:length(t_1)
            dd(3,:) = (v_2(:,2*j-1) - v_1(:,j))/(2*por_met-1);
            dd(2,:) = (dd(1,:) + dd(3,:))/2;
            for i = 1:3
                if max(abs(dd(i,:)))>m
                    m = max(abs(dd(i,:)));
                    d = dd(i,:);
                end
            end
            dd(1,:) = dd(3,:);
        end
        
        %вывод в файл количетва узлов и максимальной погрешности
        fprintf(file, '%1u %2.5g \n', length(t_2),m);
        v_1 = v_2; t_1 = t_2;
    end
    y = v_2;
    t = t_2;
end
