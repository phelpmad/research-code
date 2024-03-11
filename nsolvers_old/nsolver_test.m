L = 1; theta = 0.1; cl = 0.5; gamma = 0.1; h = 0; k = 0; flag = 0;

%[y, dy, xtrue] = mypf(L, theta, cl, gamma, h, k, flag);
[y, dy] = myBad();
xtrue = fzero(y,0);

maxit = 30; atol = 1e-8; res = 'A'; beta = 1;
m = [2 3];

a = -2; b = 2;
xl = a; xr = b;

T = mySolverTest (y, dy, xtrue, a, b, xl, xr, atol, res, maxit, beta, m, flag)

%[x, iter] = mySecant(y, a, b, atol, res, maxit, flag)



















