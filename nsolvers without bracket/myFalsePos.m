function [x, iter] = myFalsePos(f, x1, x2, atol, res, maxit, flag)

    %% choose residual type
    if res == 0 || strcmp('A', res) % absolute
        r = @(x2, x1) x1 - x2;
    elseif res == 1 || strcmp('R', res) % relative
        r = @(x2, x1) (x1 - x2)./x2;
    elseif res == 2  || strcmp('F', res) % function evaluation residual
        r = @(x2, x1) f(x2);
    end
    
    %% initialize values
    iter = 0; resnorm = 0;
    fl(1) = f(x1); fh(1) = f(x2);
    
    %% check if there is a root in the interval [x1, x2]   
    if fl*fh > 0
        fprintf("There is no root in this interval.\n");
    end

    %% set xl, xh
    if fl < 0,  xl = x1; xh = x2;
    else
        xl = x2; xh = x1;
        %% swap fl and fh values
        swap = fl;
        fl = fh;
        fh = swap;
    end
        
    %% compute interval length
    dx = xh-xl;
    
    %% start false position loop
    while 1
        rtf = xl + dx * fl/(fl-fh);
        fnew = f(rtf);
        if fnew < 0
            del = xl-rtf;
            xl = rtf;
            fl = fnew;
            resnorm = abs(del);
            x = xl;
        else
            del = xh-rtf;
            xh = rtf;
            fh = fnew;
            resnorm = abs(del);
            x = xh;
        end
        
        dx = xh-xl;
        
        if resnorm < atol || fnew == 0
            if flag ~= 0
                    fprintf('You converged in %g iterations to x = %.10f.\n', ...
                        iter, x);
            end
            break;
        elseif iter > maxit
            x = NaN; iter = iter;
            if flag ~= 0
                    fprintf('No convergence with the False Position method.\n');
            end
            break;
        else
            iter = iter+1;
        end
        
    end
 
end