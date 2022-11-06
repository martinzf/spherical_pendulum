clear; clf; clc
% Parameters
% Initial conditions
t0 = 0;
tend = 10;
theta0 = pi/2;
thetadot0 = 5;
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

% Kinetic energy
T = @(thetadot) 1/2*m*(l*thetadot).^2;
% Effective potential
U = @(theta) m*l^2*(c^2./(2*sin(theta).^2)+k*cos(theta));
% Mechanical energy
E = T(thetadot0)+U(theta0);
% Derivative of time wrt polar coordinate
dtdtheta = @(theta, t, sign) sign*1./sqrt(2*(E-U(theta)));

% Points at which U = E (turning points)
f = @(theta) U(theta)-E;
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
[t,theta] = boundedpotential(t0,theta0,sign,m,U,E,[mintheta,maxtheta],tol);
% Approximating theta(t)'s Fourier series
ftheta = fseries(t,theta,4);

% Corresponding values of phi
dphidt = @(t,phi) c./sin(ftheta(t)).^2;
[t,phi] = ode45(dphidt, [t0 tend], phi0);
% Fitting phi to polynomial
pphi = @(time) polyval(polyfit(t,phi,5),time);
fphi = @(t) mod(pphi(t),2*pi);

% Plotting
t = t0:tol:tend;
plot(t,ftheta(t),t,fphi(t))
xlabel('t(s)')
ylabel('rad')
legend('\theta', '\phi')