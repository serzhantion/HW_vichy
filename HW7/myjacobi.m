function tmp=myjacobi(k,n)
%функция возвращает массив полиномов Якоби p(i) в аналитическом виде,
%i=0,1,...,n,  i - степень полинома, k - верхний индекс.
pj=sym(zeros(n+1,1));
syms x;pj(1)=1;
switch n
    case 0
        tmp=pj;
    case 1
        pj(2)=(k+1)*x;
        tmp=pj;
    otherwise
        pj(2)=(k+1)*x;
        for i=2:n
            pj(i+1)=((i+k)/(i+2*k))*(2+(2*k-1)/i)*x*pj(i)-((i+k)/(i+2*k))*(1+(k-1)/i)*pj(i-1);
        end
        tmp=pj;
end