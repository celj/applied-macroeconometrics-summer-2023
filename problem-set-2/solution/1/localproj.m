%function beta_h = localproj(y,x,h,p)
function [beta_h, res,d12y_h]= localproj(y,x,h,p)
%%
% written by Gustavo Leyva
% Applied Macroeconometrics
% ITAM, Summer 2023
%
% constructing the estimation data, defined as ...
% the 12-month log difference of the raw time series
% y: cpi
% x: ner (nominal exchange rate)
% h: horizon in equation 1, p. 959
% p: lag; see, for instance, Figure 2, p. 974

logy = log(y);
logx = log(x);
T = length(logx);                     % length of the raw time series
for i=13:T-h
    d12y_h(i) = logy(i+h)-logy(i-12); % constructing the left-hand side variable "y_t+h" in equation 1, p. 959
    d12y(i)   = logy(i)  -logy(i-12); % used below; see lines 28 and 34
    d12x(i)   = logx(i)  -logx(i-12); % constructing the right-hand side variable "x_t" in equation 1, p. 959
end
t = 13;
g = [d12y_h(t:T-h)' d12y(t:T-h)' d12x(t:T-h)']; % grouping the constructed series
n = length(g);
if y ~= x
    X = [ones(n-p,1) g(p+1:n,3) g(1:n-p,2) g(1:n-p,3)];
    % constant, mu_h
    % "x_t" in equation 1, p. 959
    % lag of "y_t" in equation 1, p. 959
    % lag of "x_t" in equation 1, p. 959
else
    X = [ones(n-p,1) g(p+1:n,3) g(1:n-p,2)];
    % accomodates the case when "y_t" and "x_t" are the same variables
end
Y = g(p+1:n,1);
% "y_t+h" in equation 1, p. 959
beta_h = (X'*X)\(X'*Y);
res = Y - X*beta_h;
%beta_h = beta_h(2,:);% saving the coefficient associated to "x_t" in equation 1, p. 959
