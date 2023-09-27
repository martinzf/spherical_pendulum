clear; clf; clc

% Initial conditions
t0 = 0;
dt = 1e-3;
[tf,alpha0,d_alpha0,phi0,d_phi0,l,g,p_phi] = params();
t = 0:dt:tf;

% 1. Equilibrium condition, unstable or stable
if (alpha0 == 0 || alpha0 == pi) && d_alpha0 == 0
    alpha = alpha0*ones(1,length(t));
    phi = zeros(1,length(t));

% 2. Initial velocity in phi = 0 => simple pendulum, motion on a plane
elseif d_phi0 == 0
    alpha = simple_pend(alpha0,d_alpha0,l,g,t);
    phi = phi0*ones(1,length(t));

% 3. General problem
else
    [alpha, phi] = general(alpha0,d_alpha0,phi0,d_phi0,l,g,p_phi,t);
end

plot(t, alpha, t, phi)
%animate(framedur,dur,theta,phi,l)