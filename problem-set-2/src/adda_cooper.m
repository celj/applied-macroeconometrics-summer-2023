function [y_grid,p] = adda_cooper(rho,sig_u,ngrid)

sig_y = sqrt(sig_u^2/(1-rho^2));

areabreaks = linspace(0,1,ngrid);

z = norminv(areabreaks);

y_grid = z' * sig_y;

p = ones(ngrid,ngrid) / (ngrid);

% still in progress

end