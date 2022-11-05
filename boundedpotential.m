function [t,x]=boundedpotential(t0,x0,sign,m,U,E,tp,tol)
% Function solves problem of bounded motion in 1d potential
% Returns one period of motion
% t0 = initial time
% x0 = initial position
% sign = sense of integration, +-1, ie initial velocity
% m = mass
% U(x) = effective potential
% E = total energy
% tp = turning points [min max]
% tol = tolerance
    % Equation of motion
    dtdx = @(x,t) sqrt(m./(2*(E-U(x))));
    xspan = [tp(1)+tol tp(2)-tol];
    % Half motion
    [x,t] = ode45(dtdx,xspan,t0,odeset('AbsTol',tol));
    % Applying initial conditions
    x = x-x(end)+tp(2); 
    [~,idx]= min(abs(x-x0));
    t = t-sign*t(idx);
    % Full motion
    t = [t; t(end)+cumsum(flipud(diff(t)))];
    x = [x; flipud(x(1:end-1))];
end