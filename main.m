clear; clf; clc
% Parameters
% Initial conditions
t0 = 0;
tend = 60;
theta0 = pi/2;
thetadot0 = 0;
phi0 = 0;
phidot0 = 1;
% Mass
m = 1;
% Pendulum length
l = 1;
% Gravitational acceleration
g = 9.8;
% Error tolerance
tol = 1e-3;
maxiter = 100;

k = g/l;
% It can be shown L/ml^2 = sin^2(theta)phidot = const.
c = sin(theta0)^2*phidot0;

% Equilibrium condition, unstable or stable
if (theta0 < tol || pi-theta0 < tol) && abs(thetadot0) < tol
    t = t0:tol:tend;
    theta = theta0*ones(length(t),1);
    phi = zeros(length(t),1);
% Initial velocity in phi = 0 => simple pendulum, motion on a plane
elseif abs(phidot0) < tol
    % Position and velocity in theta
    thetadot = @(t,theta) [theta(2);k*sin(theta(1))];
    [t,thetavec] = ode89(thetadot, ...
        [t0 tend], ...
        [theta0;thetadot0], ...
        odeset('AbsTol',tol));
    theta = thetavec(:,1);
    % Phi
    phi = phi0*ones(length(t),1);
% General problem
else
    % Effective kinetic energy
    T = @(thetadot) 1/2*thetadot.^2;
    % Effective potential energy
    U = @(theta) c^2./(2*sin(theta).^2)+k*cos(theta);
    % Effective mechanical energy
    E = (T(thetadot0)+U(theta0))/(m*l^2);
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
        [t,theta] = boundedmotion(t0,theta0,sign, ...
            1,U,E,[mintheta,maxtheta],tol);
        % Approximating theta(t)'s Fourier series
        ftheta = fseries(t,theta,4);
    end
    dphidt = @(t,phi) c./sin(ftheta(t)).^2;
    [t,phi] = rungekutta(dphidt,t0,tend,tol, phi0);
    theta = ftheta(t);
end

% Time frames
delta = .1/tol;
t = t(1:delta:end);
theta = theta(1:delta:end);
phi = phi(1:delta:end);
% Cartesian coords
xyz = l*[sin(theta).*cos(phi);
    sin(theta).*sin(phi);
    cos(theta)];
% Line object
hold on
an = animatedline('Marker','.');
title('Spherical pendulum')
% Spherical grid
[X,Y,Z] = sphere;
X = l*X;
Y = l*Y;
Z = l*Z;
mesh(X,Y,Z)
alpha .5
hold off
% Force 3D view
view(3) 
% Axis limits
lim = 1.2*l; 
xlim("manual")
ylim("manual")
zlim("manual")
xlim([-lim lim])
ylim([-lim lim])
zlim([-lim lim])
xlabel('X')
ylabel('Y')
zlabel('Z')
grid on
% Time text
txt = annotation('textbox',[.1,.875,.1,.1],'String','t=0');
% Timer
a = tic;
for i = 1:length(t)
    clearpoints(an)
    addpoints(an,[0 xyz(1,i)],[0 xyz(2,i)],[0 xyz(3,i)])
    set(txt,'String',['t=',num2str(t(i))])
    while toc(a) < t(i)
    end
    drawnow 
end