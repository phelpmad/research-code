function [y, dy, x] = mypf(L, theta, cl, gamma, h, k, flag)
%% --------------------------------------------------------------------- %%
%% Permafrost BVP20 Semismooth function construction:
%
%   L       Latent heat
%   theta   Temperature
%   cl      Latent heat capacity
%   gamma   Controls decay of exp() on LHS
%   h       Controls hortizontal shift of function
%   k       Controls vertical shift of function;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% handle nargin
    if nargin < 7, flag = 0; end

    %% ensure correct parameter limits
    if L < 0
        L = abs(L);
    elseif L == 0
        error("L cannot equal zero!")
    end
    if cl < 0
        error("Choose cl >=0.");
    end
    if theta > 2 || theta <= 0
        error("Choose theta within (0,2].\n");
    end
    
        %% if so, make function
        y  = @(x) funct(x, L, theta, cl, gamma, h, k);
        dy = @(x) dfunct(x, L, theta, cl, gamma, h);
        x  = x0(L, theta, cl, gamma, h, k);
        
        %% plot function with flag
        if flag == 1
            b = 1; % buffer for plotting
            n = 1000; % number of points to plot with
            if x-b < 0
                t = linspace(x-b, b, n);
                plot(t,y(t),'k');
            else
                t = linspace(-b, x+b, n);
                plot(t,y(t),'k');
            end
        end
end


%% define permafrost test function
function y = funct(x,L,theta,cl,gamma,h,k)
    xval = x+h;
    y = (L*exp(theta.*xval) + gamma.*xval + k).*(xval<=0) + ...
            (cl.*xval + L + k).*(xval>0);
end

%% compute derivative of permafrost function 
function dy = dfunct(x,L,theta,cl,gamma,h)
    xval = x+h;
    dy = (theta*L.*exp(theta.*xval) + gamma).*(xval<=0) + ...
        cl.*(xval>0);
end

%% compute the true root
function x = x0(L,theta,cl,gamma,h,k)
    if k > -L
        x = fzero(@(t) L*exp(theta*t+theta*h)+gamma*t + gamma*h + k,0);
    elseif k <= -L
        x = -(k+L)/cl+h;
    end
end