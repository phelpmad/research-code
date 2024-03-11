function [x0, iter, anew, bnew] = myBisectold(f, xl, xr, atol, maxit, flag)
%% --------------------------------------------------------------------- %%
%% Bisection method: myBisect()
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
%% Tested:
%   11.10.23
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Inputs: 
%       f           function for which we want to find the root
%       xl          left endpoint of starting interval
%       xr          right endpoint of starting interval
%       atol        user specified error tolerence, i.e., stopping criteria
%       maxit       user specified maximum number of iterations
%       flag        1 prints convergence for tests
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Outputs:
%       x0      the root of f
%       iter    the number of bisections it took to get desired level
%               of accuracy for x0.
%       anew    returns new left endpoint for smaller interval
%       bnew    returns new right endpoint for smaller interval
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Example:
% f = myCf(); [x0, iter, a, b] = myBisect(f,-1,1,1e-8,10,1); 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% set initial interval [a,b] with a<b
    if xl > xr, a = xr; b = xl; else, a = xl; b = xr; end
    
    %% check if there is a root in desired interval
    if (f(a) > 0 && f(b) < 0) || (f(a) < 0 && f(b) > 0)
        if flag == 1
            fprintf('There is a root in the interval: [%g, %g]\n', a,b)
        end
    else
        error('No root in interval. Choose different values.')
    end
    
    %% bisect initial interval
    x0 = (a+b)/2; 
    
    %% initialize convergence variables
    resnorm = 0; iter = 0;
    
    %% Bisection method:
    while 1
        resnorm = abs(f(x0));    
        if resnorm < atol
            iter = iter;
            x0 = (a+b)/2; anew = a; bnew = b;
            if flag ~= 0
                fprintf('Convergence in %g iterations to x0 = %.10f\n', ...
                    iter, x0);
                fprintf('in the interval [%.10f, %.10f]\n', anew, bnew);
            end
            break;
        elseif iter > maxit
            iter = maxit;
            x0 = NaN; anew = a; bnew = b;
            if flag ~= 0
                fprintf('No convergence in %g iterations.\n', iter);
                fprintf('Try again with new values: [%.10f, %.10f]\n', ...
                    anew, bnew);
            end
            break;
        else
            if f(a) < 0
                if f(x0) < 0
                    a = x0;
                else
                    b = x0;
                end
            elseif f(a) > 0
                if f(x0) > 0
                    a = x0;
                else
                    b = x0;
                end
            end
            x0 = (a + b)/2;
            iter = iter + 1;
        end 
    end
    
end
