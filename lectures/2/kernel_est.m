%
%  KERNEL_EST.SRC
%
%  Purpose : calculate a univariate kernel (density) estimate
%            and its derivative
%  Format  : { f,d,h } = UKERNEL(xs,x,h,w,&kf);
%  Input   : xs   : (nsx1) vector, points where density is to be estimated
%            x    : (nx1) vector, observed data points
%            h    : (kx1) vector, bandwidths
%            w    : (nx1) vector, weights for different observations
%                   --> Can be used for binned data!!!
%            &kf  : pointer to procedure containing kernel function
%  Output  : f    : (nsxk) matrix, estimated density in xs, each column in f
%                   corresponds to a bandwidth
%            d    : (nsxk) matrix, derivative of estimated density
%            h    : (kx1) vector, bandwidths
%  Globals :
%  Remarks : If h<=0, the bandwidth is determined by the procedure bandws(x)
%  Written by Dick van Dijk in Gauss
%  First attempt on  7 October  1996
%  Last revision on 17 December 1997
%

function [f] = kernel_est(xs,x)

h = band_silverman(x);

% Initialization
ns = length(xs);
f = zeros(ns,1);

% Compute density estimates
i = 0;
while (i ~= ns);
    i = i+1;
    arg = (x-xs(i))./h;
    kff = k_gauss(arg);
    f(i,:) = (sum(kff)./(ns*h));
end;

end
