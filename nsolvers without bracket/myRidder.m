function [x, iter] = myRidders(f, a, b, atol, res, maxit, flag)
%% ---------------------------------------------------------------------%%
%% Ridders' method: myRidders
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Inputs:
%       f           function with desired root
%       a
%       b
%       atol        desired stopping tolerance for convergence
%       res         residual type:
%                       <'A'> absolute residual
%                       <'R'> relative residual
%                       <'F'> function residual
%       maxit       maximum number of iterations for method
%       flag        turns on or off print statements
%                       <0> turns off
%                       <1> turns on
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Outputs:
%       x               
%       iter            number of iterations to converge given tolerance
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Examples:
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% check first condition of Ridders' method
    if f(a)*f(b) > 0
        x = NaN; iter = 0;
        fprintf('There is no root within the interval [a,b]. Try again.\n')
    end

    %% choose residual type and define function r
    if strcmp('A', res) 
        r = @(x1,x2) x1 - x2;
    elseif strcmp('R', res)
        r = @(x1,x2) (x1 - x2)./x2;
    elseif strcmp('F', res)
        r = @(x1, x2) norm(f(x1)); 
    else
        error(['error:\n res = A computes absolute residual\n' ...
            ' res = R computes relative residual\n' ...
            ' res = F computes function residual\n']);
    end

    %% initialize parameters
    fl = f(a); fh = f(b);
    
    if (fl > 0 && fh < 0) || (fl < 0 && fh > 0)
        xl = a; xh = b;
        ans = 

    %% start Ridders' method
    for j = 1:maxit
        xm = 0.5*(xl+xh);
        fm = f(xm);
        s = sqrt(fm*fm-fl*fh);
        if s == 0
            x = NaN;
            fprintf('s=0 here.\n');
            break;
        end
        xnew = xm + (xm-xl)*sign(fl-fh)*fm/s;



end