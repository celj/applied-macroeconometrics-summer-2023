%
%  BANDWS.g
%
%  Purpose : calculate the optimal bandwidth in kernel estimation of a density
%            by using the rule of thumb from Silverman (1986)
%  Format  : h = BANDWS(x);
%  Input   : x   : (nx1) vector, series for which density is to be estimated
%  Output  : h   : scalar, bandwidth
%  Globals :
%  Remarks : The optimal bandwidth is calculated according to
%            eq. 3.31 of Silverman (1986)
%  Written by Dick van Dijk in Gauss
%  First attempt on  7 October  1996
%  Last revision on 15 February 1997
%
function [band] = band_silverman(x)

n = length(x);
xs = sort(x);
q1 = round(0.25*n);
q3 = round(0.75*n);
iqr = xs(q3)-xs(q1);
s = std(xs);
band = 0.9*min([s;(iqr/1.349)])/n^0.2;

end

