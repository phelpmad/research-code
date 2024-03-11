function [y, Jac, x0] = mySSf(theta,a,m,c)
%% --------------------------------------------------------------------- %%
%% Semismooth function
%
%   This needs to be further tested.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Inputs: 
%       <theta> 
%           0 << theta < 2
%           controls curve of expontential part
%           (theta = 0.1 optimal for permafrost)
%       <a>
%               controls horizontal shift
%       <m>
%           controls steepness of linear portion
%       <c>
%           c < 0 to obtain a calculated root
%           controls vertical shift
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    y = @(x) f(x,theta,m) + c; y = @(x) y(x-a);
    Jac = @(x) df(x,theta,m); Jac = @(x) Jac(x-a);
    x0 = myzero(theta,a,m,c);
end

function y = f(x,theta,m)
    y(x<=0) = exp(theta*x(x<=0));
    y(x>0) = m.*x(x>0)+1;
end

function dy = df(x,theta,m)
    dy(x<0) = theta*exp(theta*x(x<0));
    dy(x>=0) = m; % might be issues here
end

function x0 = myzero(theta,a,m,c)
    if c >= 0
        x0 = NaN;
    elseif c >= -1 && c < 0
        x0 = log(abs(c))/theta + a;
    elseif c < -1
        x0 = -(1+c)/m + a;
    end
end
