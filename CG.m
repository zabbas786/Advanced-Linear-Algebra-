function [ x, niters ] = CG( A, b, x0 )
%function [x, niters] = conjugate_gradient(A, b, x0)

% Initialization
x = x0;
r = b - A * x;
p = r;
k = 0;
niters = [];

% Iteration
while norm(r) ~= 0
    if k == 0
        p = r;
    else
        gamma = -(p' * A * r) / (p' * A * p);
        p = r + gamma * p;
    end
    alpha = (r' * r) / (p' * A * p);
    x = x + alpha * p;
    r_old = r;
    r = r - alpha * A * p;
    beta = (r' * r) / (r_old' * r_old);
    niters = [niters norm(r)];
    k = k + 1;
end

end