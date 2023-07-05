function z = simul_proj(y,x,h,p)
%%
% written by Gustavo Leyva
% Applied Macroeconometrics
% ITAM, Summer 2023
%

logy = log(y);
logx = log(x);
T = length(logx);                % length of the raw time series
for i=13:T-h
    d12y(i)   = logy(i)  -logy(i-12);
    d12x(i)   = logx(i)  -logx(i-12);
end
t = 13;
g = [d12y(t:T-h)' d12x(t:T-h)']; % grouping the constructed series
n = length(g);

if y ~= x
    X = [ones(n-p,1) g(p+1:n,2) g(1:n-p,1) g(1:n-p,2)];
    [proj_coef,proj_res]= localproj(y,x,h,p);
    proj_res_b = simple_bootstrap(proj_res);
    z = X*proj_coef + proj_res_b;
else
    X = [ones(n-p,1) g(p+1:n,2) g(1:n-p,1)];
    [proj_coef,proj_res]= localproj(y,x,h,p);
    proj_res_b = simple_bootstrap(proj_res);
    z = X*proj_coef + proj_res_b;
end
end