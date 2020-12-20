function lambda = ScalarMet(A, epsilon)
    s = size(A);
    n = s(1);
    x = rand(n, 1);
    
    lambda = (x' * A * x)/(x' * x);
    
    d = norm(A*x - lambda * x)/norm(x);
    
    i = 1;
    
    while d > epsilon 
        x = (A * x) / norm(A * x);
        d = abs(lambda-(x' * A * x)/(x' * x));
        lambda = (x' * A * x)/(x' * x);
        
        i = i + 1;
    end
     disp ('Количество итераций:');
     disp( i);
end