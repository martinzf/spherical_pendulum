function [tf,alpha0,d_alpha0,phi0,p_phi,l,g] = params()
% Function requests initial conditions
    tf = input('Duration (s): ');
    alpha0 = input('Initial south inclination (rad): ');
    d_alpha0 = input('Initial polar (up) velocity (rad/s): ');
    % If at equilibrium
    if (alpha0 == 0 || alpha0 == pi) && d_alpha0 == 0
        phi0 = 0;
        disp('Initial azimuthal coordinate = 0 rad')
    else
        phi0 = input('Initial azimuth (rad): ');
    end
    % Keep the angles in their appropriate intervals
    alpha0 = mod(alpha0,2*pi);
    if alpha0 > pi
        alpha0 = 2*pi-alpha0;
        phi0 = phi0+pi;
    end
    phi0 = mod(phi0,2*pi);
    % If bob at either pole initially, there cannot be azimuthal velocity
    if alpha0 == 0 || alpha0 == pi
        d_phi0 = 0;
        disp('Initial azimuthal velocity = 0 rad/s')
    else
        d_phi0 = input('Initial azimuthal velocity (rad/s): ');
    end
    % Pendulum length
    l = input('Pendulum length (m): ');
    % Gravitational acceleration
    g = input('Gravity (m/s^2): ');
    % Lz/ml^2 = sin^2(theta)phidot = const., Lz angular momentum wrt. Z-ax
    p_phi = sin(alpha0)^2*d_phi0;
end