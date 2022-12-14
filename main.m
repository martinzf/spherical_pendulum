clear; clf; clc

% Error tolerance
tol = 1e-3;
framedur = .1;

% Initial conditions
t0 = input('Initial time (s): ');
tf = input('Final time (s): ');
theta0 = input('Initial polar coordinate (rad): ');
thetadot0 = input('Initial polar velocity (rad/s): ');
phi0 = input('Initial azimuthal coordinate (rad): ');
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
    disp('Initial azimuthal velocity = 0 rad/s')
else
    phidot0 = input('Initial azimuthal velocity (rad/s): ');
end

% Parameters
% Pendulum length
l = input('Pendulum length (m): ');
% Gravitational acceleration
g = 9.8;
% Renamed constant
k = g/l;
% It can be shown Lz/ml^2 = sin^2(theta)phidot = const.
% where Lz is the angular momentum in the Z-axis
c = sin(theta0)^2*phidot0;

% 1. Equilibrium condition, unstable or stable
if (theta0 < tol || pi-theta0 < tol) && abs(thetadot0) < tol
    ftheta = @(t) theta0*ones(1,length(t));
    fphi = @(t) zeros(1,length(t));

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
        % Motion of theta
        [t,theta] = boundedmotion(t0,theta0,sign(thetadot0), ...
            1,U,E,[0,2*pi],tol); 
        % Interpolation to numerically obtain periodic theta(t) function
        [t,ia,~] = unique(t,'stable'); % Identify non duplicate values 
        theta = theta(ia);
        T = t(end)-t(1); % Period
        ftheta = @(tq) interp1(t,theta,t0+mod(tq-t0,T),'spline','extrap');

    % 2.2 If E = max(U) => motion in theta is bounded and non periodic
    elseif abs(E-k) < tol
        % Sense of integration
        sgn = sign(thetadot0);
        if abs(theta0) < tol && sgn == -1
            theta0 = 2*pi;
        end
        % Motion of theta
        [t,theta] = boundedmotion(t0,theta0,sgn, ...
            1,U,E,[0,2*pi],tol);
       % Interpolation for smoother animation
       [t,ia,~] = unique(t,'stable'); % Identify non duplicate values
       theta = theta(ia);
       tinf = t(end); % Time at which theta = 0 rad
       ftheta = @(tq) [interp1(t,theta,t0+tq(tq<tinf),'spline','extrap'), ...
           zeros(1,length(tq(tq>=tinf)))];

    % 2.3 If E < max(U) => bounded, periodic motion
    else
        % Find 2 points at which U = E (turning points)
        f = @(theta) U(theta)-E;
        % Find root of f using the bisection method
        mintheta = bisec(f,tol,pi,tol);
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
        t = [t;t(end)+cumsum(flipud(diff(t)))];
        theta = [theta;flipud(theta(1:end-1))];
        % Interpolation to numerically obtain periodic theta(t) function
        [t,ia,~] = unique(t,'stable'); % Identify non duplicate values 
        theta = theta(ia);
        T = t(end)-t(1); % Period
        ftheta = @(tq) interp1(t,theta,t0+mod(tq-t0,T),'spline','extrap');
    end

    % Motion of phi
    fphi = @(t) phi0*ones(1,length(t));

% 3. General problem
else
    % Effective kinetic energy
    T = @(thetadot) 1/2*thetadot.^2;
    % Effective potential energy
    U = @(theta) c^2./(2*sin(theta).^2)+k*cos(theta);
    % Effective total energy
    E = T(thetadot0)+U(theta0);
    % Slope of U
    Uprime = @(theta) -c^2*cos(theta)./sin(theta).^3-k*sin(theta);
    % Find U'(theta) = 0 => minimum 
    eqtheta = bisec(Uprime,pi/2,pi-tol,tol);

    % 3.1 Stable equilibrium in theta
    if abs(eqtheta-theta0) < tol && abs(thetadot0) < tol
        ftheta = @(t) eqtheta;

    % 3.2 Periodic motion in theta
    else
        % Find points at which U = E (turning points)
        f = @(theta) U(theta)-E;
        % Find 2 roots of f using the bisection method
        mintheta = bisec(f,tol,eqtheta,tol);
        maxtheta = bisec(f,eqtheta,pi-tol,tol);
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
        t = [t;t(end)+cumsum(flipud(diff(t)))];
        theta = [theta;flipud(theta(1:end-1))];
        % Interpolation to numerically obtain periodic theta(t) function
        [t,ia,~] = unique(t,'stable'); % Identify non duplicate values 
        theta = theta(ia);
        T = t(end)-t(1); % Period
        ftheta = @(tq) interp1(t,theta,t0+mod(tq-t0,T),'spline','extrap');
    end

    % Motion of phi
    dphidt = @(t,phi) c./sin(ftheta(t)).^2;
    [t,phi] = ode89(dphidt,[t0 tf],phi0);
    % Interpolation to numerically obtain phi(t) function
    [t,ia,~] = unique(t,'stable'); % Identify non duplicate values
    phi = phi(ia);
    T = t(end)-t(1); % Period
    fphi = @(tq) interp1(t,phi,t0+mod(tq-t0,T),'spline','extrap');
end

% Time frames
t = t0:framedur:tf;
theta = ftheta(t);
phi = fphi(t);
% Cartesian coords
xyz = l*[sin(theta).*cos(phi);
    sin(theta).*sin(phi);
    cos(theta)];
% pendulum and wake line objects
hold on
pend = animatedline('MaximumNumPoints',2,'Marker','.','MarkerSize',12, ...
    'Color','k');
wake = animatedline('MaximumNumPoints',5,'LineWidth',2, ...
    'Color',[.8 .8 .8 .5]);
title('Spherical pendulum')
% Spherical grid
[X,Y,Z] = sphere;
X = l*X;
Y = l*Y;
Z = l*Z;
mesh(X,Y,Z,'FaceAlpha','0')
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
xlabel('X (m)')
ylabel('Y (m)')
zlabel('Z (m)')
grid on
% Time text
txt = annotation('textbox',[.1,.875,.1,.1],'String','t=0');
% Timer
a = tic;
for i = 1:length(t)
    addpoints(pend,[0 xyz(1,i)],[0 xyz(2,i)],[0 xyz(3,i)])
    addpoints(wake,xyz(1,i),xyz(2,i),xyz(3,i))
    set(txt,'String',['t = ',num2str(t(i)),'s'])
    while toc(a)+t0 < t(i)
    end
    drawnow 
end