function [t,x]=boundedmotion(t0,x0,sign,m,U,E,xspan,tol)
% Function solves problem of bounded motion in 1d potential
% t0 = initial time
% x0 = initial position
% sign = sense of integration, +-1, ie initial velocity
% m = effective mass
% U(x) = effective potential
% E = effective energy
% xspan = turning points or equilibria [min max]
% tol = tolerance
    % Equation of motion
    dtdx = @(x,t) real(sqrt(m./(2*(E-U(x)))));
    xspan = [xspan(1)+tol xspan(2)-tol];
    % Half motion
    [x,t] = rungekutta(dtdx,xspan(1),xspan(2),tol,t0);
    % Applying initial conditions 
    [~,idx]= min(abs(x-x0));
    t = t-sign*t(idx)+t0;
end