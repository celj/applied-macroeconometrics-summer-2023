function c = sole_spline(x,fx)

A = [1 x(1) 0 0 0;1 x(2) 0 0 0;0 0 1 x(2) x(2)^2;0 0 1 x(3) x(3)^2;0 1 0 -1 -2*x(2)];
b = [fx(1);fx(2);fx(2);fx(3);0];
c = gaussj(A,b);

end