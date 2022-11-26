clear; clf; clc

% Initial conditions
t0 = 0;
tend = 20;
theta0 = pi/2;
thetadot0 = 4.427;
phi0 = 0;
phidot0 = 0;

% Parameters
% Mass
m = 1;
% Pendulum length
l = 1;
% Gravitational acceleration
g = 9.8;
% Renamed constant
k = g/l;
% It can be shown L/ml^2 = sin^2(theta)phidot = const.
c = sin(theta0)^2*phidot0;
% Error tolerance
tol = 1e-3;
maxiter = 100;

% Keep the angles in their appropriate intervals
theta0 = mod(theta0,2*pi);
if theta0 > pi
    theta0 = 2*pi-theta0;
    phi0 = phi0+pi;
end
phi0 = mod(phi0,2*pi);
% If the mass is at either of the poles initially, it cannot have azimuthal
% velocity
if abs(theta0) < tol || abs(theta0-pi) < tol
    phidot0 = 0;
end

% 1. Equilibrium condition, unstable or stable
if (theta0 < tol || pi-theta0 < tol) && abs(thetadot0) < tol
    t = t0:tol:tend;
    theta = theta0*ones(length(t),1);
    phi = zeros(length(t),1);

% 2. Initial velocity in phi = 0 => simple pendulum, motion on a plane
elseif abs(phidot0) < tol
    % Effective kinetic energy
    T = @(thetadot) 1/2*thetadot.^2;
    % Effective potential
    U = @(theta) k*cos(theta);
    % Effective total energy
    E = T(thetadot0)+U(theta0);

    % 2.1 If E > max(U) => motion in theta is unbounded
    if E > k+tol
        thetadot = @(t,theta) sign(thetadot0)*sqrt(2*(E-U(theta)));
        [t,theta] = rungekutta(thetadot,t0,tend,tol,theta0);

    % 2.2 If E = max(U) => motion in theta is bounded and non periodic
    elseif abs(E-k) < tol
        sgn = sign(thetadot0);
        if abs(theta0) < tol && sgn == -1
            theta0 = 2*pi;
        end
        [t,theta] = boundedmotion(t0,theta0,sgn, ...
            1,U,E,[0,2*pi],tol);
        t = [t t(end)+tol:tol:tend];
        theta = [theta zeros(1,length(t)-length(theta))];

    % 2.3 If E < max(U) => bounded, periodic motion
    else
        % Find 2 points at which U = E (turning points)
        f = @(theta) U(theta)-E;
        % Find root of f using the bisection method
        mintheta = bisec(f,tol,pi,tol,maxiter);
        maxtheta = 2*pi-mintheta;
        % Which turning point are we moving towards?
        % Non-negative velocity
        if thetadot0 ~= 0
            sgn = sign(thetadot0);
            % Located at left turning point
        elseif abs(theta0-mintheta) < abs(theta0-maxtheta)
            sgn = 1;
            % Located at right turning point
        else
            sgn = -1;
        end
        % Motion of theta
        [t,theta] = boundedmotion(t0,theta0,sgn, ...
            1,U,E,[mintheta,maxtheta],tol); % One period
        % Full period
        t = [t,t(end)+cumsum(fliplr(diff(t)))];
        theta = [theta,fliplr(theta(1,1:end-1))];
        % Approximating theta(t)'s Fourier series
        ftheta = fseries(t,theta,4);
        t = t0:tol:tend;
        theta = ftheta(t);
    end

    % Motion of phi
    phi = phi0*ones(1,length(t));

% 3. General problem
else
    % Effective kinetic energy
    T = @(thetadot) 1/2*thetadot.^2;
    % Effective potential energy
    U = @(theta) c^2./(2*sin(theta).^2)+k*cos(theta);
    % Effective total energy
    E = (T(thetadot0)+U(theta0))/(m*l^2);
    % Slope of U
    Uprime = @(theta) -c^2*cos(theta)./sin(theta).^3-k*sin(theta);
    % Find U'(theta) = 0 => minimum 
    eqtheta = bisec(Uprime,pi/2,pi-tol,tol,maxiter);

    % 3.1 Stable equilibrium in theta
    if abs(eqtheta-theta0) < tol && abs(thetadot0) < tol
        ftheta = @(t) eqtheta;

    % 3.2 Periodic motion in theta
    else
        % Find points at which U = E (turning points)
        f = @(theta) U(theta)-E;
        % Find 2 roots of f using the bisection method
        mintheta = bisec(f,tol,eqtheta,tol,maxiter);
        maxtheta = bisec(f,eqtheta,pi-tol,tol,maxiter);
        % Which turning point are we moving towards?
        % Non-negative velocity
        if thetadot0 ~= 0
            sgn = sgn(thetadot0);
            % Located at left turning point
        elseif abs(theta0-mintheta) < abs(theta0-maxtheta)
            sgn = 1;
            % Located at right turning point
        else
            sgn = -1;
        end
        % Motion of theta
        [t,theta] = boundedmotion(t0,theta0,sgn, ...
            1,U,E,[mintheta,maxtheta],tol); % One period
        % Full period
        t = [t,t(end)+cumsum(fliplr(diff(t)))];
        theta = [theta,fliplr(theta(1,1:end-1))];
        % Approximating theta(t)'s Fourier series
        ftheta = fseries(t,theta,4);
    end

    % Motion of phi
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
alpha .1
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