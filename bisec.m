function c = bisec(f,a,b,tol,maxiter)
% Calculates the root of a function in an interval [a,b] by bisection
    it=1;
    % Bolzano's theorem => need f(a) to have a different sign from f(b)
    if sign(f(a))==sign(f(b))
        warning('f(a) and f(b) have the same sign')
        c=NaN;
    else
        % Midpoint between a and b
        c=(a+b)/2;
        % Check sign of f(c) and repeat procedure with c the other point
        % whose sign is different from c's 

        % While the algorithm hasn't converged and we haven't reached
        % maximum iterations
        while abs(f(c)) > tol && it < maxiter
            it=it+1;
            if sign(f(c))==sign(f(a))
                a=c;
            else % if sign(fc) == sign(fb)
                b=c;
            end
            c=(a+b)/2;
        end
        % Non convergence warning
        if abs(f(c))>tol && it==maxiter
            warning('Did not converge with given max iterations')
            disp(['Answer given with an absolute error of ',num2str(abs(f(c)))])
        end
    end
end