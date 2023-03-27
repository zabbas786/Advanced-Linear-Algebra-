function [x, niters] = Method_of_Steepest_Descent(A, b, x0)
    r = b - A * x0;
    p = r;
    niters = 0;
    while norm(r) ~= 0
        q = A * p;
        alpha = dot(p, r) / dot(p, q);
        x = x0 + alpha * p;
        r = r - alpha * q;
        p = r;
        x0 = x;
        niters = niters + 1;
    end
end