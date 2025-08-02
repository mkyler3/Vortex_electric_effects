function [x, solR, solW] = electric_effects(l, gamma, m, xmax)
h = xmax / m; % initiate grid
N = m - 2;
xmin = h / 2;
x = (xmin:h:xmax)';

b = 5;  % initial condition

solW = [0; -b; zeros(N,1)];
while min(solW) < solW(1) % since exact b.c. unknown, ensure minumum of W occurs at x=0
    b = min(solW);
    guessW = [-b * ones(N+1, 1); 0]; % initial guesses, including boundary conditions
    guessR = [0; ones(N+1, 1)];
    solW = guessW;
    solR = guessR;
    W1 = guessW(2:end-1); % remove endpoints to keep boundaries fixed at boundary conditions 
    R1 = guessR(2:end-1);
    full_sol = [W1; R1];

    g = 1 ./ x; % inhomogeneous part of the equations
    g1 = g(2:end-1);

    J_TL = sparse(N,N); % set up Jacobian pieces, construct full matrix later
    J_TR = sparse(N,N);
    J_BL = sparse(N,N);
    J_BR = sparse(N,N);
    
    V_T = zeros(N, 1); % set up vector that multiplies inverse Jacobian
    V_B = zeros(N, 1); 

    tolerance = 1e-10; % set tolerance, error, and max iterations of Newton's method
    error = 10;
    iter = 0;
    iter_max = 1000;

    while error > tolerance && iter < iter_max % Newton's method

        for i=1:N % Construct Jacobian and V vector components
            J_TR(i,i) = (2)*(1/(gamma^2))*R1(i);
            J_BL(i,i) = R1(i); 
            J_BR(i,i) = (l*g1(i))^2 + (2/(h^2)) + 3*(R1(i)^2) - 1 + W1(i);
            V_B(i) = (l*g1(i))^2 * solR(i+1) - (g1(i)/(2*h))*(solR(i+2)-solR(i))...
                - (1/(h^2)) * (solR(i+2)-2*solR(i+1) + solR(i))...
                + (solR(i+1))^3 - solR(i+1) + solW(i+1)*solR(i+1);
        if (i < N)
            J_TL(i,i+1) = 1/(h^2) + (g1(i))/(2*h);  
            J_BR(i,i+1) = (- (g1(i))/(2*h) - (1/(h^2)));
        end
        if (i > 1)
            V_T(i) = (g1(i)/(2*h))*(solW(i+2)-solW(i)) + (1/(h^2))*(solW(i+2)...
                -2*solW(i+1)+solW(i)) + (1/gamma^2)*((solR(i+1))^2 - 1);
            J_TL(i,i) = - (2/(h^2)) ;
            J_TL(i,i-1) = - (g1(i)/(2*h)) + (1/(h^2));
            J_BR(i,i-1) = ((g1(i)/(2*h)) - (1/h^2));
        end
        if i == 1
            V_T(i) = (g1(i)/2*h)*(solW(i+2)-solW(i+1)+(((g(1)/(2*h))...
                + 1/(h^2))^-1)*(1/gamma^2)) + (1/h^2)*(solW(i+2)-solW(i+1)-(((g(1)/(2*h))...
                + 1/(h^2))^-1 )*(1/gamma^2)) + (1/gamma^2)*(solR(i+1)^2 -1);
            J_TL(i,i) = -g1(i)/(2*h) - 1/(h^2);
        end
    
        end
        J = sparse([J_TL J_TR; J_BL J_BR]); % construct full Jacobian
        V = [V_T; V_B];

        delta = -J \ V; % calculate change to solution
        full_sol = full_sol + delta; % update full solution
        error = max(abs(delta)); % update error
        iter = iter + 1;

        new_w1 = full_sol(1) - 1/((g(1)/(2*h)) + 1/(h^2)); % enforce symmetry condition at zero
        solW = [new_w1; full_sol(1:N); 0]; 
        solR = [0; full_sol(N+1:end); 1];
        W1 = full_sol(1:N);
        R1 = full_sol(N+1:end);
    end
end
end