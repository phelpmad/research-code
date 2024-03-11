function [x, iter] = myNewtonAnderson1(y, dy, beta, guess, atol, res, maxit, flag)
%% ---------------------------------------------------------------------%%
%% Newton-Anderson (1)
%   
%   This algorithm is the Newton-Anderson (1) method for scalar functions
%   using Algorithm 1 in Pollock RAM 2020 paper. This algorithm initializes
%   with Newtons method and then performs an update step called the
%   Anderson-step in the following code using a dampening parameter beta.
%   For the tests in the paper, this value should be equal to 1 and more
%   tests should be performed with beta ~= 1.
%
%   The value w keeps track of Newton's update step and gamma is the
%   Anderson parameter computed at each step.
%
%% Tested: 
%       11.10.23        used myCf() and beta = 1 but has not been tested
%                       with beta not equal to 1.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Inputs:
%       y           function with desired root
%       dy          Jacobian of function
%       beta        uniform dampening parameter
%       atol        desired stopping tolerance for convergence
%       res         residual type:
%                       <'A'> absolute residual
%                       <'R'> relative residual
%       maxit       maximum number of iterations for method
%       flag        turns on or off print statements
%                       <0> turns off
%                       <1> turns on
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Outputs:
%       x           returns the root of the function f = y
%       iter        returns the total number of iterations used to converge
%                   or maxit if reached.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Examples:
% [y, dy] = myCf();
% [x, iter] = myNewtonAnderson1(y, dy, 1, 1, 1e-8, 'A', 10, 1);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% choose residual type and define function r
    if strcmp('A', res) 
        r = @(x1,x2) x1 - x2;
    elseif strcmp('R', res)
        r = @(x1,x2) (x1 - x2)./x2;
    elseif strcmp('F', res)
        r = @(x1,x2) norm(x1);
    else
        error(['error:\n res = A computes absolute residual\n' ...
            ' res = R computes relative residual\n']);
    end

    %% initialize variables
    iter = 0; resnorm = 0; 
    xold = guess;
    w = - dy(xold) \ y(xold);
    x = xold + beta * w;
    
    %% start Newton-Anderson (1) method
    while 1

        wnew = - dy(x) \ y(x);                         % Newton's step
        gamma = (wnew' *(wnew - w))/norm(wnew - w).^2; % Anderson step
        xnew = x + beta*wnew ...
                - gamma * ( (x - xold) + beta * (wnew - w) );
        resnorm = abs(r(xnew,x));
    
        if resnorm < atol 
            iter = iter; xold = x; x = xnew;
            if flag ~=0 
                fprintf("Convergence in %g iterations to x0 = %.12f\n", iter, x);
            end
            break;
        elseif iter > maxit
            iter = iter;
            x = NaN;
            if flag ~= 0
                fprintf("No convergence. The maximum number of iterations were exceeded.\n");
            end
            break;
        else
            iter = iter+1; xold = x; x = xnew; w = wnew;
        end
        %
    end
end
