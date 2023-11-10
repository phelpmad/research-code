function [iter, soln] = newtonanderson_m (x0, m, beta)
%% Newton-Anderson(m)
%
%% ---------------------------------------------------------------------%%
%% Parameters:
% <x0>                                % choose initial guess
% <m>                                 % sets AA window depth
% <beta> = 1;                         % uniform damping parameter in (0,1]
%% ---------------------------------------------------------------------%%
%
    %% set up default case for scalar AA with m = 2
    if nargin < 1, x0 = 1; m = 2; beta = 1;
    elseif nargin < 2, m = 2; beta = 1; 
    elseif nargin < 3, beta = 1; 
    end

    %% set tolerances
    atol = 1e-8; maxit = 100; 
    
    %% set up function and its jacobian
    [f, Jac] = myfunct();
 
    %% initialize variables
    iter = 0; resnorm = 0;
    %x = zeros(maxit+1,1); w = zeros(maxit+1,1);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Initialize newton anderson m starting values
    w(1) = - Jac(x0) \ f(x0);
    x(1) = x0 + beta * w(1);
    %%
    %% Start for loop k = 1,2,3,...
    for k = 1:maxit

        %% Compute next update step w
        w(k+1) = - Jac(x(k)) \ f(x(k));

        %% Set anderson window size mk
        mk = min(k,m);
        
        %% Take differences of update steps and iterates
        %% Note: 
            % The last values of W and X should be disregarded, and
            % they are not used in following calculations.
        W = flip(w(1:end)) - flip([0, w(1:end-1)]);
        X = flip(x(1:end)) - flip([0, x(1:end-1)]);
        
        %% Compute F and E
        % chooses the first mk entries of W and X
        F = W(1:mk); E = X(1:mk);

        %% Compute gamma
        gamma = F \ w(k+1);
        
        %% uncomment to see mk, w, F, gamma
        % mk, w(k+1), F, gamma, pause

        %% Compute next iterate
        x(k+1) = x(k) + beta * w(k+1) - (E + beta .* F) * gamma;

        %% Compute residual 
        resnorm = abs(f(x(k+1)));    
  
        %% Check convergence
        if resnorm < atol
            iter = iter;
            soln = x(k+1);
            fprintf("converged in %g iterations to x* = %.12f\n", iter, x(k+1));
            break;
        elseif iter > maxit
            iter = iter; soln = 0;
            break;
        else
            iter = iter+1;
        end
    end
%%
end
%
%% ---------------------------------------------------------------------%%
%% functions:
function [y, dy] = myfunct()
    y = @(x) exp(-x) - x;
    dy = @(x) -exp(-x) - 1;
end

%% calculates argmin as in Alg. 3
function v = argmin(n,w,F)
    r = rand(10000,n);
    if n == 1
        %% scalar case
        for i = 1:length(r), rnorm(i) = norm(r(i),2); end
        u = r(find(rnorm<=1));
        for i = 1:length(u), x(i) = norm(w - F*u(i),2); end
        [~,idx] = min(x);
        v = u(idx);
    else 
        %% for n > 1
        % calculate 2-norm of randomly generated vectors
        for i = 1:length(r), rnorm(i) = norm(r(i,:),2); end
        % find vectors with magnitude leq 1
        u = r(find(rnorm<=1),:);
        % compute w-Fu
        for i = 1:length(u), x(i) = norm(w - F*u(i,:),2); end
        % find index that minimizes norm
        [~,idx] = min(x);
        % return argmin
        v = u(idx,:);
    end
end
%% ---------------------------------------------------------------------%%
