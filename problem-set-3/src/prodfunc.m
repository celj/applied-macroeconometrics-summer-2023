function [F] = prodfunc(k)

global A alpha delta

F = A * k .^ alpha + (1 - delta) .* k;

end