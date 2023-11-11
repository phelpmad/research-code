%% Contiuous function used to test nsolvers
function [y, Jac] = myCf()
    y = @(x) exp(-x)-x;
    Jac = @(x) -exp(-x)-1;
end