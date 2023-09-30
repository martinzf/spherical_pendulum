clear; clc

% Initial conditions
t0 = 0;
dt = 1/30;
[tf,alpha0,d_alpha0,phi0,p_phi,l,g] = params();
t = 0:dt:tf;

% 1. Equilibrium condition, unstable or stable
if (alpha0 == 0 || alpha0 == pi) && d_alpha0 == 0
    alpha = alpha0*ones(1,length(t));
    phi = zeros(1,length(t));

% 2. Initial velocity in phi = 0 => simple pendulum, motion on a plane
elseif p_phi == 0
    alpha = simple_pend(alpha0,d_alpha0,l,g,t);
    phi = phi0*ones(1,length(t));

% 3. General problem
else
    [alpha, phi] = general(alpha0,d_alpha0,phi0,p_phi,l,g,t);
end

animate(dt,tf,alpha,phi,l)