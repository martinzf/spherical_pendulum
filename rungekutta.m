function [x,y]=rungekutta(dydx,a,b,dx,y0)
    %Runge-Kutta approximation of ODE's
    y = y0; 
    x = a;
    it = 1; 
    while x(it)<b
        y = [y y(:,it)+ ...
            dx*dydx(x(it)+dx/2,y(:,it)+ ...
            dx/2*dydx(x(it),y(:,it)))];
        x(it+1) = x(it)+dx;
        it = it+1;
    end
end