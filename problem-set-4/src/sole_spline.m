function [c] = sole_spline(x,fx)

c = zeros(5,1);
n=size(c,1);
normVal=Inf;
tol=1e-5;

A = [
    1 0 x(2) 0 0;
    0 1 0 x(3) (x(3)^2);
    0 0 1 0 0;
    0 0 -((x(3) + x(2))/(x(3) - x(2))) 1 0;
    0 0 0 (1/(x(2) + x(3))) 1;
    ];

b = [
    fx(2);
    fx(3);
    (fx(1) - fx(2))/(x(1) - x(2));
    (2 * x(2) * (fx(2) - fx(3)) / ((x(2) - x(3))^2));
    ((fx(2) - fx(3)) / ((x(2) - x(3))^2));
    ];

% Gauss Seidel Method
while normVal>tol
    c_old=c;
    for i=1:n
        sigma=0;
        for j=1:i-1
            sigma=sigma+A(i,j)*c(j);
        end
        for j=i+1:n
            sigma=sigma+A(i,j)*c_old(j);
        end
        c(i)=(1/A(i,i))*(b(i)-sigma);
    end
    normVal=norm(c_old-c);
end

end