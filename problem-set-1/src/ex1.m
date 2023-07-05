tol = 1e-6;
kmax = 1e4;

funList = {@(x) cos(x) - x, @(x) atan(x), @(x) sin(x), @(x) x.^2, @(x) x.^2 + 1};
funListLatex = {'$\cos(x) - x$', '$\arctan(x)$', '$\sin(x)$', '$x^2$', '$x^2 + 1$'};

x0 = [0.5, 1, 3, 0.5, 10];

for i = 1:length(funList)
    fprintf('Function %d\n', i)
    fprintf('--------------------\n')
    
    [iterations, solutions, errors] = newton(funList{i}, x0(i), tol, kmax);
    
    figure;
    plot(iterations, errors);
    title(sprintf('Function %d: %s', i, funListLatex{i}));
    subtitle(sprintf('Newton method, $x_0 = %.4f$, $e_T = %.4f$', x0(i), errors(end)));
    xlabel('Iterations');
    ylabel('Error');
    
    figure;
    plot(iterations, solutions);
    title(sprintf('Function %d: %s', i, funListLatex{i}));
    subtitle(sprintf('Newton method, $x_0 = %.4f$, $x_T = %.4f$', x0(i), solutions(end)));
    xlabel('Iterations');
    ylabel('Solution');
    
    fprintf('Newton method\n Iterations = %d,\n Solution = %.4f,\n Error = %.4f\n', iterations(end), solutions(end), errors(end))
    
    [iterations, solutions, errors] = chord(funList{i}, x0(i), tol, kmax);
    
    figure;
    plot(iterations, errors);
    title(sprintf('Function %d: %s', i, funListLatex{i}));
    subtitle(sprintf('Chord method, $x_0 = %.4f$, $e_T = %.4f$', x0(i), errors(end)));
    xlabel('Iterations');
    ylabel('Error');
    
    figure;
    plot(iterations, solutions);
    title(sprintf('Function %d: %s', i, funListLatex{i}));
    subtitle(sprintf('Chord method, $x_0 = %.4f$, $x_T = %.4f$', x0(i), solutions(end)));
    xlabel('Iterations');
    ylabel('Solution');
    
    fprintf('Chord method\n Iterations = %d,\n Solution = %.4f,\n Error = %.4f\n', iterations(end), solutions(end), errors(end))
    
    fprintf('\n')
end
