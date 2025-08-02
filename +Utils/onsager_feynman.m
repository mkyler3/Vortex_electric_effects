function [x, sol] = onsager_feynman(l, m, xmax)
h = xmax / m; % initiate grid
N = m - 2;
xmin = h / 2;
x = (xmin:h:xmax)';
g = 1 ./ x;
g1 = g(2:end-1);

sol = [0; ones(N+1, 1)];
sol1 = sol(2:end-1);
J = sparse(N, N);
V = zeros(N, 1);

tolerance = 1e-10; % set tolerance, error, and max iterations of Newton's method
error = 10; 
iter = 0;
iter_max = 1000;

while error > tolerance && iter < iter_max % Newton's method
    for i = 1:N % construct Jacobian and V vector components
        J(i, i) = g1(i) * (1/h) - (2 / h^2) - (l*g1(i))^2 - 3 * sol1(i)^2 + 1;
        V(i) = (g1(i)/h)*(sol(i+1) - sol(i)) + (1/h^2)*(sol(i+2) - 2*sol(i+1) + sol(i)) ...
             - (l*g1(i))^2 * sol(i+1) - sol(i+1)^3 + sol(i+1);
    end
    for i = 1:N
        if i < N
            J(i, i+1) = 1 / h^2;
        end
        if i > 1
            J(i, i-1) = -g1(i)*(1/h) + 1/h^2;
        end
    end
    delta = -J \ V; % calculate change to solution
    sol1 = sol1 + delta; % update solution
    error = max(abs(delta)); % update error
    sol = [0; sol1; 1];
    iter = iter + 1;
end
end