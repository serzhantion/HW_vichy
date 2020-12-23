function [phi, dphi, ddphi]=myphi_dphi_d2kk_good_anal(k,n)
%функция возвращает массивы коодинатных функций phi(i), i=1,...,n
%в аналитическом виде и массивы их производных.
phi=sym(zeros(n,1));
dphi=sym(zeros(n,1));
ddphi=sym(zeros(n,1));
syms x;
jac=myjacobi(k,n-1);
%массив многочленов Якоби p(i), i=0,...,n-1
%в аналитическом виде,  k - верхний индекс.
djac=myjacobi(k-1,n);
for i=1:n
    phi(i)=(1-x^2)*jac(i);
    phi(i)=collect(phi(i));
    %приведение многочлена к виду, содержащему степени
    %с соответствующими коэффициентами.
    dphi(i)=-2*(i)*(1-x^2)^(k-1)*djac(i+1);
    dphi(i)=collect(dphi(i));
    ddphi(i)=-2*(i)*((k-1)*(1-x^2)^(k-2)*(-2*x)*djac(i+1)+(1-x^2)^(k-1)*(i+2*(k-1)+1)/2)*jac(i);
    ddphi(i)=collect(ddphi(i));
end
end