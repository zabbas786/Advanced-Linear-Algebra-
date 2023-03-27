function [x, niters] = PCG(A, b, x0)
% Initialization
x = x0;
r = b - A * x;
M = diag(diag(A));
Minv = inv(M);
z = Minv * r;
p = z;
niters = 0;

% Iteration
while norm(r) ~= 0
    if niters == 0
        p = z;
    else
        gamma = -((p' * A * r) / (p' * A * p));
        p = z + gamma * p;
    end
    alpha = (r' * z) / (p' * A * p);
    x = x + alpha * p;
    r = r - alpha * A * p;
    z = Minv * r;
    niters = niters + 1;
end