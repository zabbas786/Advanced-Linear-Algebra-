function [ x, niters ] = Method_of_Steepest_Descent_ichol( A, b, x0 )

A = ichol(sparse(A), struct('type','ict','droptol',1e-3,'michol','off'));

% Method of Steepest Descent for solving Ax = b
% A: matrix
% b: column vector
% x0: initial guess for x
% tol: tolerance for stopping criteria
% maxIter: maximum number of iterations allowed
% x: solution vector
% k: number of iterations used

r = b - A * x0;
k = 0;
x = x0;
niters = 0;
while norm(r) ~= 0
    p = r;
    q = A * p;
    alpha = (p' * r) / (p' * q);
    x = x + alpha * p;
    r = r - alpha * q;
    niters = niters + 1;
    k = k + 1;
end
end