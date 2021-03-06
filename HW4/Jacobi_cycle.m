function lambda = Jacobi_cycle(A, eps)
    A_abs_v_triu = abs(triu(A,1));
    siz = size(A);
    n = siz(1);
    
    num = 0;
    m = max(max(A_abs_v_triu));
    while m > eps
        for I_m = 1:n
            for J_m = I_m+1:n
            V = eye(siz);
            phi = atan(-2*A(I_m,J_m)/(A(J_m,J_m) - A(I_m,I_m)))/2;
            c = cos(phi);
            s = sin(phi); 
        
            V(I_m,J_m) = -s;
            V(J_m,I_m) = s;
               
            V(I_m,I_m) = c;
            V(J_m,J_m) = c;
        
            A =A*V;
        
            A=transpose(V)*A;     
            A_abs_v_triu = abs(triu(A,1));
            m = max(max(A_abs_v_triu));
        
            num = num + 1;
            end
        end
    end
    disp ('Количество итераций:');
    disp(num);
    lambda = diag(A);
end
