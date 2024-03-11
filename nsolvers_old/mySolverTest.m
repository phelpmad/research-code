function T = mySolverTest (y, dy, xtrue, a, b, xl, xr, atol, res, maxit, beta, m, flag) 
    
    %% choose guesses
    guess = linspace(xtrue-abs(xtrue)/2, xtrue+abs(xtrue)/2, 5);
    
    if nargin < 11
        beta = 1; m = 2; flag = 0;
    end
    
    %% create names for table rows
    name = {'Newton', 'Secant', 'Bisect', 'Bisect-Newton', ...
                'Newton-Anderson 1'};
    for i = 1:length(m)
        name(end+1) = {sprintf('Newton-Anderson %g',m(i))};
    end
    
    for i = 1:length(guess)
        [x(:,i), iter(:,i)] = generateSolverTest(y, dy, guess(i), a, b, xl, xr, atol, res, maxit, beta, m, flag) ;
    end

    %% generate table based on guess to the left and right of xtrue root
    k = 1; n = 2*length(guess);
    for j = 1:2:n
        A(:,j) = x(:,k); varname(j) = {sprintf('x, guess = %g', guess(k))};
        A(:,j+1) = iter(:,k); varname(j+1) = {sprintf('iter%g', k)};
        k = k+1;
    end
    
    T = array2table(A, 'RowNames', name, 'VariableNames', varname);

end

function [x, iter] = generateSolverTest(y, dy, guess, a, b, xl, xr, atol, res, maxit, beta, m, flag) 
    
    [x(1), iter(1)]   = myNewton(y, dy, guess, atol, res, maxit, flag);
    [x(2), iter(2)]   = mySecant(y, a, b, atol, res, maxit, flag);
    [x(3), iter(3)]   = myBisect(y, xl, xr, atol, maxit, flag);
    [x(4), iter(4)]   = myBisectNewton(y, dy, xl, xr, atol, res, maxit, flag);
    [x(5), iter(5)]   = myNewtonAnderson1(y, dy, beta, guess, atol, res, maxit, flag);

    for i = 1:length(m)
        [x(5+i), iter(5+i)] = myNewtonAndersonm(y, dy, m(i), beta, guess, atol, res, maxit, flag);
    end

end
