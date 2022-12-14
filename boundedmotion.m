function [t,x]=boundedmotion(t0,x0,sign,m,U,E,xspan,tol)
% Function solves problem of bounded motion in 1d potential
% t0 = initial time
% x0 = initial position
% sign = sense of integration, +-1, ie initial velocity
% m = effective mass
% U(x) = effective potential
% E = effective energy
% xspan = bounds of motion
% tol = tolerance
    % Equation of motion
    dtdx = @(x,t) real(sqrt(m./(2*(E-U(x)))));
    xspan = [xspan(1)+tol xspan(2)-tol];
    % Half motion
    [x,t] = ode89(dtdx,xspan,0);
    % Applying initial conditions 
    [~,idx]= min(abs(x-x0));
    t = t-sign*t(idx)+t0;
end