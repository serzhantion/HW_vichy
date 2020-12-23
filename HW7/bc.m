function [q_r, q_l] = bc(alpha,beta, a, b, p, y, z)
fun = p*y*z;
if or(alpha(1)==0,alpha(2)==0)
    q_l = 0;
else
    q_l = double(subs(fun,'x',a))*alpha(1)/alpha(2);
end
if or(beta(1)==0,beta(2)==0)
    q_r = 0;
else
    q_r = double(subs(fun,'x',b))*beta(1)/beta(2);
end
end