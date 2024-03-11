function [x, iter] = myBroyden(f, A, guess, atol, res, maxit, flag)
%% ---------------------------------------------------------------------%%
%% Broyden's Method: 
% 
%   Multidimenisonal generalization of the Secant method in higher
%   dimensions. Start with initial guess x0 = [x1; x2; ... ; xn] and use A
%   as the scaled identity (tbd) and iterate with:
%       x_k+1 = x_k - A_k \ f(x_k)
%       A_k+1 = A_k + f(x_k+1)(x_k+1-x_k)^T / norm(x_k+1-x_k)^2
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Inputs:
%       f           function
%       guess       initial starting point
%       atol        absolute error tolerance
%       maxit       maximum iterations desired
%       res         error function
%       <res>
%                   0 || 'A' - absolute error: abs(b - a)
%                   1 || 'R' - relative error: abs((a - b)/b)
%                   2 - residual error: abs(f(b))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Examples:
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% set argument presets
    if nargin < 6, flag = 0; end
    
    %% choose residual type
    if res == 0 || strcmp('A', res) % absolute
        r = @(x2, x1) x1 - x2;
    elseif res == 1 || strcmp('R', res) % relative
        r = @(x2, x1) (x1 - x2)./x2;
    elseif res == 2 || strcmp('F', res) % residual
        r = @(x2, x1) norm(f(x2));
    end

    %% Initialize values
    iter = 0; resnorm = 0;
    x = guess;    

    %% Secant method
    while 1
        %% computer steps
        xnew = x - A \ f(x);
        Anew = A + f(xnew)*(xnew-x)' / ((xnew-x)'*(xnew-x));

        %% compute residual
        resnorm = abs(r(xnew, x));
        
        %% check convergence
        if resnorm < atol
            iter = iter; x = xnew;
            if flag ~= 0
                fprintf('You converged in %g iterations with Broydens method to\n x = [', ...
                    iter);
                fprintf('%.10f ', x);
                fprintf(']\n')
            end
            break;
        elseif iter > maxit
            iter = iter; x = NaN;
            if flag ~= 0
                fprintf('No convergence with Broydens method.\n');
            end
            break;
        else
            x = xnew; A = Anew;
            iter = iter + 1;
        end
    end
  
end
