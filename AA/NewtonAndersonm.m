%% Newton-Anderson(m)
% Verified with m = 2 and F,E hard coded
%% ---------------------------------------------------------------------%%
%% Parameters:
m = 2;                              % choose depth
beta = 1;                           % uniform damping parameter in (0,1]
guess = 1;                          % choose initial guess
[f, Jac] = myfunct();               % set function with its derivative
atol = 1e-8; maxit = 10;           % set stopping criteria
iter = 0; resnorm = 0;              % initialize counters
%% ---------------------------------------------------------------------%%
%
x(1) = guess; 
w(2) = - Jac(x(1)) \ f(x(1));
x(2) = x(1) + beta * w(2);
%
for k = 2:maxit
    w(k+1) = - Jac(x(k)) \ f(x(k)), pause
    mk = min(k-1,m);
    W = w(k-1:k+1);
    if mk < m
        F = W(2)-W(1);
        E = x(k)-x(k-1);
    else
        X = x(k-2:k);
        F = [W(3)-W(2), W(2)-W(1)];
        E = [X(3)-X(2), X(2)-X(1)];
    end
    %gamma = inv(F'*F)*F'*w(k+2)
    %gamma = lsqr(F,w(k+1));
    gamma = F\w(k+1), pause
    x(k+1) = x(k) + beta*w(k+1) ...
            - (E + beta*F)*gamma;

    % compute residual based on Pollock RAM paper
    %resnorm = abs(f(xnew));    
    resnorm = abs(x(k+1) - fzero(f,guess));
    %resnorm = abs(x(k+1) - x(k));

    if resnorm < atol % converged
        iter = iter; %xold = x; x = xnew;
        fprintf("converged in %g iterations to x0 = %.12f\n", k, x(end));
        break;
    elseif iter > maxit % exceeded max iterations
        break;
    else
        iter = iter+1; %xold = x; x = xnew;
    end

end
%
%% ---------------------------------------------------------------------%%
%% functions:
function [y, dy] = myfunct()
    y = @(x) exp(-x) - x;
    dy = @(x) -exp(-x) - 1;
end
%% ---------------------------------------------------------------------%%
