function c=bisec(f,a,b,tol,maxiter)
    it=1;
    if sign(f(a))==sign(f(b))
        warning('f(a) and f(b) have the same sign')
        c=NaN;
    else
        c=(a+b)/2;
        while abs(f(c))>tol && it<maxiter
            it=it+1;
            if sign(f(c))==sign(f(a))
                a=c;
            else %if sign(fc)==sign(fb)
                b=c;
            end
            c=(a+b)/2;
        end
        if abs(f(c))>tol && it==maxiter
            warning('Did not converge with given max iterations')
            disp(['Answer given with an absolute error of ',num2str(abs(f(c)))])
        end
    end
end