function [x, iter] = myNewton(y, dy, guess, atol, res, maxit, flag)
%% ---------------------------------------------------------------------%%
%% Newtons Method: myNewton()
%
%   Use Newton's method to compute the root of a given scalar function.
%   This function requires the function and its derivative, a guess, a
%   tolerance for convergence, a maximum number of iterations and the flag
%   command will turn on (0) or off (1) print statements.
%
%% Tested: 11.10.23 with smooth continuous function myCf()
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Inputs:
%       y           function with desired root
%       dy          Jacobian of function
%       guess       initial iterate for Newton's method
%       atol        desired stopping tolerance for convergence
%       res         residual type:
%                       <'A'> absolute residual
%                       <'R'> relative residual
%       maxit       maximum number of iterations for method
%       flag        turns on or off print statements
%                       <0> turns off
%                       <1> turns on
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Outputs:
%       x               Newton's method solution
%       iter            number of iterations to converge given tolerance
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Examples:
% [y, dy] = myCf(); 
% guess = 0.1; atol = 1e-8; res = 'A'; maxit = 10; flag = 1;
% [x, iter] = myNewton(y, dy, guess, atol, res, maxit, flag);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% set dampening constant (hard coded)
    beta = 1;

    %% choose residual type and define function r
    if strcmp('A', res) 
        r = @(x1,x2) x1 - x2;
    elseif strcmp('R', res)
        r = @(x1,x2) (x1 - x2)./x2;
    elseif strcmp('F', res)
        r = @(x1, x2) norm(y(x1)); 
    else
        error(['error:\n res = A computes absolute residual\n' ...
            ' res = R computes relative residual\n']);
    end

    %% initialize parameters
    x = guess; iter = 0; resnorm = 0;
    
    %% start Newton's Method
    while 1
        Jac = dy(x);
        f = y(x);
        xnew = x - beta * Jac \ f;
        resnorm = abs(r(xnew,x));
    
        if resnorm < atol
            iter = iter;
            x = xnew;
            if flag ~= 0
                fprintf('You converged in %g iterations with Newtons method to\n x = [', ...
                    iter);
                fprintf('%.10f ', x);
                fprintf(']\n')
            end
            break;
        elseif iter > maxit
            iter = iter; x = NaN;
            if flag ~= 0
                fprintf('No convergence.\n Iterations: %g\n', iter-1);
            end
            break;
        else
            iter = iter + 1;
            x = xnew;
        end
    end
end