function [t,x]=boundedmotion(t0,x0,sign,m,U,E,tp,tol)
% Function solves problem of bounded motion in 1d potential
% Returns one period of motion
% t0 = initial time
% x0 = initial position
% sign = sense of integration, +-1, ie initial velocity
% m = effective mass
% U(x) = effective potential
% E = effective energy
% tp = turning points [min max]
% tol = tolerance
    % Equation of motion
    dtdx = @(x,t) sqrt(m./(2*(E-U(x))));
    xspan = [tp(1)+tol tp(2)-tol];
    % Half motion
    [x,t] = ode89(dtdx,xspan,t0,odeset('AbsTol',tol));
    % Applying initial conditions 
    [~,idx]= min(abs(x-x0));
    t = t-sign*t(idx)+t0;
    % Full motion
    t = [t; t(end)+cumsum(flipud(diff(t)))];
    x = [x; flipud(x(1:end-1))];
end