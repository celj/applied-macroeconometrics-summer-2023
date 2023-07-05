function [iter_store,sol_store,funct_store] = chord(x0)
%%
% CHORD.M
% Applied Macroeconometrics
% ITAM, Summer 2023
% Written by Gustavo Leyva
%

global ii
crit = 100;
tol = 10^-5;
i = 0;
sol_store(1) = x0;
iter_store(1) = i;
eval(['funct_store(1) = norm(func' num2str(ii) '(x0),inf);']);
%funct_store(1) = norm(func3(x0),'inf');
eval(['[f0,df0] = func' num2str(ii) '(x0);']);
%[f0,df0] = func3(x0);
while crit > tol
    i = i+1;
    iter_store(i+1) = i;
    eval(['[f,df] = func' num2str(ii) '(x0);']);
    %[f,df] = func3(x0);
    dirtn = -df0\f;
    x1 = x0+dirtn;
    crit = norm(x1-x0);
    sol_store(i+1) = x1;
    eval(['funct_store(i+1) = norm(func' num2str(ii) '(x1),inf)./funct_store(1);']);
    %funct_store(i+1) = norm(func3(x1),'inf')./funct_store(1);
    x0 = x1;
end