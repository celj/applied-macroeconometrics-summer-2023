function[betas_store] = bootstrap(y,x,h_max,p,betas,simulations)

betas_store = zeros(simulations,length(betas));
u = y - x * betas;

for i = 1:simulations
    y_resid = x * betas + simple_bootstrap(u);
    for j = 0:h_max
        betas_store(i,j+1) = localproj(y_resid,x,j,p);
    end
end

end