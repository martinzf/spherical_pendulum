function [alpha,phi] = general(alpha0,d_alpha0,phi0,p_phi,l,g,t)
    % Equilibrium in alpha
    if sin(alpha0)^4/cos(alpha0) == p_phi^2*l/g
        alpha = alpha0*ones(length(t));
        phi = mod(phi0 + p_phi/sin(alpha0)^2*t, 2*pi);
    % Generic case
    else
        % Energy integral
        h = d_alpha0^2/2 + p_phi^2/(2*sin(alpha0)^2) - g/l*cos(alpha0);
        % change of variables polynomial
        p = [-g/l, -h, g/l, h - p_phi^2/2];
        e = sort(roots(p)); % Roots in ascending order
        e1 = e(3);
        e2 = e(2);
        e3 = e(1);
        % Sign for all +- equations
        sgn = initial_direction(d_alpha0,alpha0,acos(e2),acos(e1));
        % Solution for alpha
        lambda = sqrt(e1-e3)/2;
        m1 = (e1-e2)/(e1-e3);
        e0 = 0;
        if sgn < 0
            e0 = ellipticK(m1);
        end
        if sign(d_alpha0)
            e0 = 2*ellipticK(m1) + ...
                ellipticF(asin(sqrt(abs((cos(alpha0)-e2)/((cos(alpha0)-e3)*m1)))), m1);
        end
        [sn,~,~] = ellipj(-sgn*lambda*t*sqrt(2*g/l)+e0,m1); % Expression evaluated at (t-t0)
        alpha = acos((m1*e3*sn.^2-e2)./(m1*sn.^2-1));
        % Solution for phi
        phi = phi0 + cumtrapz(t,p_phi./sin(alpha).^2);
        phi = mod(phi, 2*pi);
    end
end