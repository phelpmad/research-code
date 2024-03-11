%[y, iter] = myNewtonAndersonm(f, Jac, m, beta, x0, atol, res, maxit, flag)
L = 1; theta = 0.1; cl = 0.5; gamma = 0.1; h = 0; k = 0; flag = 0;

[y, dy, xtrue] = mypf(L, theta, cl, gamma, h, k, flag);

atol = 1e-8; res = 'A'; beta = 1;