clear; clf; clc
% Parameters
% Initial conditions
t0 = 0;
tend = 1;
theta0 = pi;
thetadot0 = 0;
phi0 = 0;
phidot0 = 2;
% Mass
m = 1;
% Gravitational acceleration
g = 9.8;
% Pendulum length
l = 1;
% Error tolerance
tol = 1e-3;
maxiter = 100;

k = g/l;
% It can be shown L/ml^2 = sin^2(theta)phidot = const.
c = sin(theta0)^2*phidot0;

% Equilibrium condition, stable or unstable
if (theta0 < tol || pi-theta0 < tol) && abs(thetadot0) < tol
    ftheta = @(t) theta0*ones(length(t),1);
    fphi = @(t) zeros(length(t),1);
% Initial velocity in phi = 0 => simple pendulum, motion on a plane
elseif abs(phidot0) < tol
    fphi = @(t) phi0*ones(length(t),1);
    % Position and velocity in theta
    thetadot = @(t,theta) [theta(2);k*sin(theta(1))];
    [t,theta] = ode45(thetadot,[t0 tend],[theta0;thetadot0]);
% General problem
else
    % Kinetic energy
    T = @(thetadot) 1/2*m*(l*thetadot).^2;
    % Effective potential
    U = @(theta) m*l^2*(c^2./(2*sin(theta).^2)+k*cos(theta));
    % Mechanical energy
    E = T(thetadot0)+U(theta0);
    % Slope of U
    Uprime = @(theta) -c^2*cos(theta)./sin(theta).^3-k*sin(theta);
    % Find U'(theta) = 0 => minimum 
    eqtheta = bisec(Uprime,pi/2,pi-tol,tol,maxiter);
    % Stable equilibrium in theta
    if abs(eqtheta-theta0) < tol && abs(thetadot0) < tol
        ftheta = @(t) eqtheta;
    else
        % Find points at which U = E (turning points)
        f = @(theta) U(theta)-E;
        % Find 2 roots of f using the bisection method
        mintheta = bisec(f,tol,eqtheta,tol,maxiter);
        maxtheta = bisec(f,eqtheta,pi-tol,tol,maxiter);
        % Which turning point are we moving towards?
        % Positive velocity
        if thetadot0 > 0
            sign = 1;
            % Negative velocity
        elseif thetadot0 < 0
            sign = -1;
            % Located at left turning point
        elseif abs(theta0-mintheta) < abs(theta0-maxtheta)
            sign = 1;
            % Located at right turning point
        else
            sign = -1;
        end
        % Motion of theta
        [t,theta] = boundedpotential(t0,theta0,sign,m,U,E,[mintheta,maxtheta],tol);
        % Approximating theta(t)'s Fourier series
        ftheta = fseries(t,theta,4);
    end
    dphidt = @(t,phi) c./sin(ftheta(t)).^2;
    [t,phi] = ode45(dphidt, [t0 tend], phi0,odeset('AbsTol',tol));
    % Fitting phi to polynomial
    pphi = @(time) polyval(polyfit(t,phi,5),time);
    fphi = @(t) mod(pphi(t),2*pi);
end

% Plotting
t = t0:tol:tend;
plot(t,ftheta(t),t,fphi(t))
xlabel('t(s)')
ylabel('rad')
legend('\theta', '\phi')