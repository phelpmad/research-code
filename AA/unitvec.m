r = rand(10000,2);
for i = 1:length(r)
    rnorm(i) = norm(r(i,:),2);
end
idr = find(rnorm<=1);
u = r(find(rnorm<=1),:);
for i = 1:length(u)
    mval(i) = 2-[0.44 -0.002]*u(i,:)';
end
[argmin, idx] = min(mval);

gamma = u(idx,:)
