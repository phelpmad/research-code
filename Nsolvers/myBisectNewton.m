function [x, iter] = myBisectNewton(y, dy, xl, xr, atol, res, maxit, flag)
%% --------------------------------------------------------------------- %%
%% Bisection method with Newton's method: myBisectNewton()
%
%   The bisection method is a root finding algorithm that takes a function
%   and an interval and first checks if there is a root in the desired
%   interval. For example, if f(a) < 0 and f(b) > 0 then f must be equal to
%   zero at some point in the interval [a,b].
%
%   This method checks if there is a root in the interval, also makes sure
%   the interval is in the correct ordering of a<b to avoid an infinite
%   loop.
%
%   At each step we check the absolute value of f(x0) to determine if the
%   value is within a given absolute tolerance <atol>. 
%
%   The interval [a,b] is bisected at each step and we check the sign of
%   f(x) at the endpoints to determine the new interval, e.g., the new 
%   iterval is either [(a+b)/2, b] or [a, (a+b)/2]. We continues till a 
%   maximum iteration is met or the method converged.
%
%   If the maximum number of iterations is reached, then we use Newton's
%   method to continue to find the root if possible where the guess is the
%   final value reached using the Bisection method. The root of f is
%   returned as a single value of x. Then, iter is the total number of
%   iterations needed to converge with tolerance <atol>. 
%
%   This code needs parameters from myBisect() and from myNewton().
%   
%   Within the code, x keeps track of the roots found from each method. By
%   first keeping track of the root of f using the Bisection method with
%   x(1) and then x(2) gives the root found using Newton's method.
%   Similarly, iter keeps track of the respective iterations used in each.
%
%% Tested: 11.10.23 with myCf()
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Inputs: 
%       y           function for which we want to find the root
%       dy          derivative of the function for Newton's method
%       xl          left endpoint of starting interval
%       xr          right endpoint of starting interval
%       atol        user specified error tolerence, i.e., stopping criteria
%       res         residual type:
%                       <'A'> absolute residual
%                       <'R'> relative residual
%       maxit       user specified maximum number of iterations
%       flag        turns on or off print statements
%                       <0> turns off
%                       <1> turns on
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Outputs:
%       x          returns the final value of the root using first the
%                  Bisection method and second (if applicable) Newton's
%                  method. 
%       iter       returns the total number of iterations used within the
%                  algorithm which equals the number of iterations from 
%                  the Bisection method and (if applicable) Newton's
%                  method.   
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Example:
% [y, dy] = myCf();
% [x, iter] = myBisectNewton(y, dy, -1, 1, 1e-8, 'A', 5, 1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% set initial interval [a,b] with a<b
    if xl > xr, a = xr; b = xl; else, a = xl; b = xr; end
    
    %% check if there is a root in desired interval
    if (y(a) > 0 && y(b) < 0) || (y(a) < 0 && y(b) > 0)
        if flag == 1
            fprintf('There is a root in the interval: [%g, %g]\n', a,b)
        end
    else
        error('No root in interval. Choose different values.')
    end
    
    %% compute midpoint of initial interval [a,b]
    x0 = (a+b)/2; 
    
    %% initialize convergence variables
    resnorm = 0; iter = [0 0]; x = [0 0];
    
    %% Bisection method:
    while 1
        resnorm = abs(y(x0));    
        if resnorm < atol
            iter = sum(iter);
            x = (a+b)/2; anew = a; bnew = b;
            if flag == 1
                fprintf('Convergence in %g iterations to x0 = %.10f\n', ...
                    iter, x);
                fprintf('in the interval [%.10f, %.10f]\n', anew, bnew);
            end
            break;
        elseif iter(1) > maxit
            iter(1) = maxit;
            x(1) = (a+b)/2;
            if flag ~=0
                fprintf('No convergence with the Bisection method.\n');
                fprintf('Using Newtons method with guess = %g\n', x(1));
            end
            [x(2), iter(2)] = myNewton(y, dy, x(1), atol, res, maxit, 0);
            if iter(2) > maxit && flag ~= 0
                fprintf('No convergence using Newtons method\n');
            elseif flag ~=0
                fprintf('Convergence in %g iterations with Newtons method.\n', iter(2));
                fprintf('Total number of iterations is %g.\n', sum(iter));
            end
            iter = sum(iter); % return total number of iterations
            x = x(2); % return final result
            break;
        else
            if y(a) < 0
                if y(x0) < 0
                    a = x0;
                else
                    b = x0;
                end
            elseif y(a) > 0
                if y(x0) > 0
                    a = x0;
                else
                    b = x0;
                end
            end
            x0 = (a + b)/2;
            iter(1) = iter(1) + 1;
        end 
    end
    
end
