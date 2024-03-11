function [y, iter] = myNewtonAndersonm(f, Jac, m, beta, guess, atol, res, maxit, flag)
%% ---------------------------------------------------------------------%%
%% Newton-Anderson (m)
%   
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Inputs:
%       f           function with desired root
%       Jac         Jacobian of function
%       m           Anderson window depth
%       beta        uniform dampening parameter
%       x0          initial guess for Newton's method set up
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
% [f, Jac] = myCf();
% m = 2; beta = 1; guess = 1; 
% atol = 1e-8; res = 'A'; maxit = 10; flag = 0;
% [x, iter] = myNewtonAndersonm(f, Jac, m, beta, guess, atol, res, maxit, flag);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
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

    %% initialize variables
    iter = 0; resnorm = 0;

%% --------------------------------------------------------------------- %%
%% Scalar Newton-Anderson (m) method
%%
    if length(guess) == 1

        w(1) = - Jac(guess) \ f(guess);
        x(1) = guess + beta * w(1);
        
        %% start Newton-Anderson (m) method
        for k = 1:100
    
            %% Compute next update step w
            w(k+1) = - Jac(x(k)) \ f(x(k));
    
            %% Set anderson window size mk
            mk = min(k,m);
            
            %% Take differences of update steps and iterates
            W = flip(w(1:end)) - flip([0, w(1:end-1)]);
            X = flip(x(1:end)) - flip([0, x(1:end-1)]);
            
            %% Compute F and E
            % chooses the first mk entries of W and X
            F = W(1:mk); E = X(1:mk);
    
            %% Compute gamma 
            gamma = F \ w(k+1);
            
            %% Compute next iterate
            x(k+1) = x(k) + beta * w(k+1) - (E + beta .* F) * gamma;
    
            %% Compute residual 
            resnorm = abs(r(x(k+1), x(k)));   
    
            while 1
                %% Check convergence
                if resnorm < atol
                    iter = iter; y = x(k+1); done = true;
                    if flag ~= 0
                        fprintf("Convergence in %g iterations to x* = %.12f\n", iter, y)
                    end
                    break;
                elseif iter > maxit
                    iter = iter; y= NaN; done = true;
                    if flag ~= 0 
                        fprintf("No convergence using Newton-Anderson(m).\n")
                    end
                    break;
                else
                iter = iter + 1; done = false;
                break; % continue in for loop
                end
            end
    
            %% exit for loop if convergence or maxit reached
            if done == true
                break;
            end
        end

%% -------------------------------------------------------------------- %%    
%% Multidimensional Newton-Anderson (m) method
%%
    else
        w(:,1) = - Jac(guess) \ f(guess);
        x(:,1) = guess + beta * w(:,1);
        
        %% start Newton-Anderson (m) method
        for k = 1:100
        
            %% Compute next update step w
            w(:,k+1) = - Jac(x(:,k)) \ f(x(:,k));
        
            %% Set anderson window size mk
            mk = min(k,m);
            
            %% Take differences of update steps and iterates
            W = flip(w(:,1:end)) - flip([zeros(length(guess),1), w(:,1:end-1)]);
            X = flip(x(:,1:end)) - flip([zeros(length(guess),1), x(:,1:end-1)]);
            
            %% Compute F and E
            % chooses the first mk entries of W and X
            F = W(:,1:mk); E = X(:,1:mk);
        
            %% Compute gamma 
            gamma = F \ w(:,k+1);
            
            %% Compute next iterate
            x(:,k+1) = x(:,k) + beta * w(:,k+1) - (E + beta .* F) * gamma;
        
            %% Compute residual 
            resnorm = norm(f(x(:,k+1)));   
        
            while 1
                %% Check convergence
                if resnorm < atol
                    iter = iter; y = x(:,k+1); done = true;
                    if flag ~= 0
                        fprintf('You converged in %g iterations with Newton-Anderson method to\n x = [', ...
                            iter);
                        fprintf('%.10f ', y);
                        fprintf(']\n')
                    end
                    break;
                elseif iter > maxit
                    iter = iter; y= x(:,k+1); done = true;
                    if flag ~= 0 
                        fprintf("No convergence using Newton-Anderson(m).\n")
                    end
                    break;
                else
                iter = iter + 1; done = false;
                break; % continue in for loop
                end
            end
        
            %% exit for loop if convergence or maxit reached
            if done == true
                break;
            end
        end
    end
%% --------------------------------------------------------------------- %%
end