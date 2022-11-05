function x = newton(f,fprime,x0,tol,maxiter)
% Root finding function which uses Newton's method
% f = function of interest (single variable)
% fprime = derivative of f
% x0 = initial guess of root
% tol = tolerance in f
% maxiter = maximum number of iterations
    it = 1;
    x = x0-f(x0)/fprime(x0);
    while abs(f(x)) > tol && it < maxiter
        x0 = x;
        x = x0-f(x0)/fprime(x0);
        it = it+1;
     end
     if abs(f(x)) > tol && it == maxiter
        warning('Did not converge with given max iterations')
        disp(['Answer given with an absolute error of ',num2str(abs(f(x)))])
     end
end