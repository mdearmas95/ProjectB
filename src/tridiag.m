function x= tridiag(e,f,g,r)
%Tridiagonal matrix solver
%input
%e= subdiagonal vector 
%   f = diagonal vector 
%   g = superdiagonal vector 
%   r = right hand side vector 
% output: 
%   x = solution vector

n = length(f);
x = r;
% forward elimination
for k = 2:n-1
    factor = e(k-1)/f(k-1);
    f(k) = f(k) - factor*g(k-1);
    r(k) = r(k) - factor*r(k-1);
end

% back substitution
x(n) = r(n)/f(n);
for k = n-1:-1:2
    x(k) = (r(k)-g(k)*x(k+1))/f(k);
end