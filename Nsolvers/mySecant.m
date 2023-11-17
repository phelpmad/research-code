function [x, iter] = mySecant(f, a, b, atol, res, maxit, flag)
%% ---------------------------------------------------------------------%%
%% Secant Method: mySecant()
% 
%   Use the secant method to find a root. This function takes in a 
%   function and two initial points for which we can define the secant
%   line to approximate the derivative as an approximation of Newton's 
%   method as a root finding algorithm.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Inputs:
%       f           function
%       a          first starting point
%       b          second starting point
%       atol        absolute error tolerance
%       maxit       maximum iterations desired
%       res         error function
%       <res>
%                   0 || 'A' - absolute error: abs(b - a)
%                   1 || 'R' - relative error: abs((a - b)/b)
%                   2 - residual error: abs(f(b))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Examples:
% f = myCf(); a = 0; b = 1; atol = 1e-8; res = 'A'; maxit = 20; flag = 1;
% [x, iter] = mySecant(f, a, b, atol, res, maxit, flag);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% set argument presets
    if nargin < 5, maxit = 20; restype = 0; end
    if nargin < 6, restype = 0; end
    
    %% choose residual type
    if res == 0 || strcmp('A', res) % absolute
        r = @(x2, x1) x1 - x2;
    elseif res == 1 || strcmp('R', res) % relative
        r = @(x2, x1) (x1 - x2)./x2;
    elseif res == 2 % residual
        r = @(x2, x1) f(x2);
    end

    %% Initialize values
    iter = 0; resnorm = 0;
    x(1) = a; x(2) = b;
    x(3) = x(2) - f(x(2))*(x(2) - x(1))./(f(x(2)) - f(x(1)));
    
    %% Secant method
    for k = 3:100
        
        %% compute next value
        x(k+1) = x(k) - f(x(k))*(x(k) - x(k-1))./(f(x(k)) - f(x(k-1)));

        %% compute residual
        resnorm = abs(r(x(k+1), x(k)));
        
        %% check convergence
        while 1
            if resnorm < atol
                iter = iter; x = x(k+1); done = true;
                if flag ~= 0
                    fprintf('You converged in %g iterations to x0 = %.10f.\n', ...
                        iter, x);
                end
                break;
            elseif iter > maxit
                iter = iter; x = NaN; done = true;
                if flag ~= 0
                    fprintf('No convergence with Secant method.\n');
                end
                break;
            else
                done = false;
                iter = iter + 1;
            end
        end

        %% exit for loop if method converged or maxit reached
        if done == true
            break;
        end

    end
  
end
