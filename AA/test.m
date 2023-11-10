clear all

n = [1 5 2 7 23];

for k = 2:4
    mk = min(k-1,2)
    N = n(k-1:k+1)
    if mk < 2
        F = N(2)-N(1)
    else
        F = [N(3) - N(2), N(2) - N(1)]
    end
end