%% Newton's Method
% 10.04.2023 -- verified this code works nicely
%% ---------------------------------------------------------------------%%
%% Parameters:
beta = 1;                           % uniform damping parameter in (0,1]
guess = 0;                          % choose initial guess
[y, dy] = myfunct();                % set function with its derivative
atol = 1e-8; maxit = 10;           % set stopping criteria
x = guess; iter = 0; resnorm = 0;   % initialize
%% ---------------------------------------------------------------------%%
%
while 1
    Jac = dy(x); % evaluate Jac at current guess
    f = y(x); % evaluate f at current guess

    % Newton's Step:
    xnew = x - beta * Jac \ f;

    % calculate residual
    resnorm = abs(x - xnew); % absolute
    %resnorm = abs((x - xnew)./xnew); % relative

    % Iterate
    if resnorm < atol % converged
        iter = iter; x = xnew;
        fprintf("converged in %g iterations to x0 = %.12f\n", iter, x);
        break;
    elseif iter > maxit % exceeded max iterations
        break;
    else
        iter = iter+1; x = xnew;
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
