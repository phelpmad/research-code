%% Newton-Anderson(1)
% VERIFIED with f(x) = exp(-x)-x
% plot these functions
%% ---------------------------------------------------------------------%%
%% Parameters:
beta = 1;                           % uniform damping parameter in (0,1]
guess = 1;                          % choose initial guess
[f, Jac] = myfunct();                % set function with its derivative
atol = 1e-8; maxit = 100;            % set stopping criteria
iter = 0; resnorm = 0;              % initialize counters
%% ---------------------------------------------------------------------%%
%
xold = guess; 
w = - Jac(xold) \ f(xold);
x = xold + beta * w;
%
while 1
    wnew = - Jac(x) \ f(x);                         % Newton's step
    gamma = abs(wnew - w) \ (wnew *(wnew - w));     % Anderson's step
    
    xnew = x + beta*wnew ...
            - gamma * ( (x - xold) + beta * (wnew - w) );

    % compute residual based on Pollock RAM paper
    %resnorm = abs(f(xnew));    
    %resnorm = abs(xnew - fzero(f,guess));
    resnorm = abs(xnew - x);

    if resnorm < atol % converged
        iter = iter; xold = x; x = xnew;
        fprintf("converged in %g iterations to x0 = %.12f\n", iter, x);
        break;
    elseif iter > maxit % exceeded max iterations
        break;
    else
        iter = iter+1; xold = x; x = xnew;
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
