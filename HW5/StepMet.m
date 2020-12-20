function lambda = StepMet(A, epsilon)
    s = size(A);
    n = s(1);
    x = rand(n, 1);
    x(1)=max(x);
    x_t = A * x;
    
    x = x_t;
    
    lambda = norm(x_t(1)) / norm(x(1));
    
    d = 1;
    
    i = 1;
    
    while d > epsilon
        d = lambda;
        x_t = A * x;
        lambda = norm(x_t(1)) / norm(x(1));
        d = abs(d-lambda);
        i = i + 1;
        x = x_t / norm(A * x);
    end
    disp ('Количество итераций:');
    disp( i);
end