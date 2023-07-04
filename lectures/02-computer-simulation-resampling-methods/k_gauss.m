function [k] = k_gauss(u)

k = pdf('norm',u,0,1);

end