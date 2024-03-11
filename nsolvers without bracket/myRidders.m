function [x, iter] = myRidders(f, x1, x2, atol, res, maxit, flag)
%% Ridders' Method:
%
%   A root finding algorithm that combines the false position method with
%   bisection method to bracket the root. Given x1 and x2 such that f has a
%   root between these values, we determine the midpoint x3 = 0.5*(x1+x2)
%   and then evaulate f at all three points. (Still not sold on this yet)
%   But Ridders' method then solves uniquely for an exponential function
%   that then turns the residual into a straight line. This accelerates the
%   root finding algorithm. The exponential function is determined using
%   the quadratic formula depending on f(xi) for i = 1,2,3. The update step
%   for xnew is guaranteed to lie within the interval, hence, bracketing
%   the solution (for example, Newton's method does not bracket the
%   solution and singularities in the derivative can cause f(xi) to tend to
%   infinity). 

%% ------------------------------------------------------------------- %%

    %% choose residual type and define function r
    if strcmp('A', res) 
        r = @(x1,x2) x1 - x2;
    elseif strcmp('R', res)
        r = @(x1,x2) (x1 - x2)./x2;
    elseif strcmp('F', res)
        r = @(x1,x2) f(x1);
    else
        error(['error:\n res = A computes absolute residual\n' ...
            ' res = R computes relative residual\n']);
    end
    
    %% set iteration counter
    iter = 0; resnorm = 0;

    %% initialize values
    fl = f(x1);
    fh = f(x2);

    if (fl > 0 && fh < 0) || (fl < 0 && fh > 0)
        xl = x1;
        xh = x2;
        x = -1.1e30; % unlikely value 

        while 1

            %% first function evaluation in each iteration
            xm = 0.5*(xl+xh);
            fm = f(xm);
            s = sqrt(fm*fm-fl*fh);

            %% exit if s = 0
            if s == 0
                if flag ~= 0
                    fprintf("Cannot divide by zero.\n")
                end
                break;
            end

            %% compute update step and resnorm
            xnew = xm + (xm-xl) * ( s \ sign(fl-fh) * fm );
            resnorm = norm(r(xnew, x));
    
            %% check convergence in first function evaulation
            if resnorm < atol
                x = xnew; iter = iter;
                if flag ~= 0
                    fprintf("You converged using Ridders method " + ...
                        "during the first function evauluation.\n");
                end
                break;
            end
            
            %% second function evaluation in each iteration
            x = xnew;
            fnew = f(x);

            %% check if we found the true root in second function evaluation
            if fnew == 0
                x = xnew; iter = iter;
                if flag ~= 0
                    fprintf("You converged using Ridders method" + ...
                        "during the second function evaluation.\n");
                end
                break;
            end

            %% bookkeeping for next iteration to keep root bracketed
            if sign(fnew) ~= sign(fm)
                xl = xm;
                fl = fm;
                xh = x;
                fh = fnew;
            elseif sign(fnew) ~= sign(fl)
                xh = x;
                fh = fnew;
            elseif sign(fnew) ~= sign(fh)
                xl = x;
                fl = fnew;
            else
                fprintf("The root is not bracketed.\n");
                break;
            end
            
            %% compute resnorm from second function evaluation
            resnorm = norm(r(xh,xl));

            %% check convergence in second function evaluation
            if resnorm < atol
                x = xnew; iter = iter;
                if flag ~= 0
                    fprintf("You converged using Ridders method" + ...
                        "during the second function evaluation.\n");
                end
                break;
            %% break if maxit reached
            elseif iter > maxit
                iter = iter; x = NaN;
                if flag ~= 0
                    fprintf("No convergence.\n")
                end
                break;
            %% continue to next iteration
            else
                iter = iter + 1;
            end
        end
    elseif fl == 0
        x = x1;
        if flag ~= 0
            fprintf("The given x1 value is the true root.\n")
        end
    elseif fh == 0
        x = x2;
        if flag ~= 0
            fprintf("The given x2 value is the true root.\n")
        end
    else
        error("The root must be bracketed.")
    end

end