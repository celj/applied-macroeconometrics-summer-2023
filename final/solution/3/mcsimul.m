function [mcs] = mcsimul(s,p,m)
%%
% MCSIMUL.M
% Numerical methods course
% Summer 2013
% Written by Gustavo Leyva
% University of Minnesota
% Department of Economics
% gustavo.leyva.jimenez@gmail.com
% Created on 07.15.13
% Modified on 07.17.13
% 
% PURPOSE   Simulates a Markov chain with states s and stochastic kernel p
% USAGE     [mcs] = mcsimul(s,p,m)
% INPUTS    s       : states (kx1)
%           p       : stochastic kernel (Markov transition matrix) (kxk)
%           m       : length of the simulation (1x1)
% OUTPUTS   mcs     : simulated series of length m and associated
%                     indexes (2x1)
% REFERENCES
%           1.Häggström, O. (2003). "Finite Markov Chains and Algorithmic
%           Applications." London Mathematical Society. Student Texts 00.
%           Cambridge University Press.
%    
ergo = ergodic(p);
mcs = zeros(m,2);

for i = 1:m;
    if i == 1;
        [x0,ind] = ini_func(rand(1,1),ergo,s);
        mcs(i,:) = [x0 ind];
    else
        [x0,ind] = upd_func(rand(1,1),p,s,x0);        
        mcs(i,:) = [x0 ind];
        x0 = mcs(i,1);
    end;    
end;

end