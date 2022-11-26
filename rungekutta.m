function [x,y] = rungekutta(dydx,a,b,dx,y0)
% 2nd order Runge-Kutta ODE numerical solver
% dydx = analytical expression for the derivative of y wrt x
% [a,b] = integration bounds
% dx = step size in x
% y0 = initial y value
    y = y0; 
    x = a;
    it = 1; 
    % While the algorithm hasn't converged
    while x(it)<b
        y = [y y(:,it)+ ...
            dx*dydx(x(it)+dx/2,y(:,it)+ ...
            dx/2*dydx(x(it),y(:,it)))];
        x(it+1) = x(it)+dx;
        it = it+1;
    end
end