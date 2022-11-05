clear; clf; clc
% Parameters
% Initial conditions
t0 = 0;
tend = 1;
theta0 = pi/2;
thetadot0 = 0;
phi0 = 0;
phidot0 = 1;
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

% Kinetic energy / mass * length^2
T = @(thetadot) 1/2*thetadot.^2;
% Effective potential
U = @(theta) c^2./(2*sin(theta).^2)+k*cos(theta);
% Mechanical energy / mass * length^2
E = T(thetadot0)+U(theta0);
% Derivative of time wrt polar coordinate
dtdtheta = @(theta, t, sign) sign*1./sqrt(2*(E-U(theta)));

% Points at which U = E (turning points)
f = @(theta) U(theta) - E;
% Derivative of f
f1 = @(theta) -c^2*cos(theta)./sin(theta).^3-k*sin(theta);
mintheta = newton(f,f1,1,tol,maxiter);
maxtheta = newton(f,f1,3,tol,maxiter);

% Equilibrium => theta = 0, theta = pi, theta = stable min
% Find stable equilibrium
% Second derivative of f
f2 = @(theta) c^2*(cos(2*theta)+2)./sin(theta).^4-k*cos(theta);
eqtheta = newton(f1,f2,3*pi/4,tol,maxiter);
% Check equilibrium condition
eqpos = theta0 < tol || pi-theta0 < tol || abs(theta0-eqtheta) < tol;
if eqpos && thetadot0 == 0
    return
end

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
[t,theta] = boundedpotential(t0,theta0,sign,m,U,E,[mintheta, maxtheta], tol);
% Fitting theta(t)
thetafit = polyfit(t,theta,4);

% Corresponding values of phi
dtdphi = @(phi,t) sin(polyval(thetafit, t)).^2/c;
[phi, t] = ode45(dtdphi, [0 2*pi], t0);
% Fitting phi
phifit = polyfit(t,phi,4);

% Plotting
t = t0:tol:tend;
plot(t,polyval(thetafit,t),t,polyval(phifit,t))
xlabel('t(s)')
ylabel('rad')
legend('\theta', '\phi')