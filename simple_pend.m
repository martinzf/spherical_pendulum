function alpha = simple_pend(alpha0,d_alpha0,l,g,t)
% Check Wikipedia for arbitrary amplitude solution of simple pendulum
% With slight manipulation one can find a solution for d_theta0 =/= 0
% It's easier to find solutions for h >= g/l
    % Energy integral
    h = d_alpha0^2/2 - g/l*cos(alpha0);
    % Sign for all +- equations
    sgn = initial_direction(d_alpha0,alpha0,-pi,pi);
    % Bounded periodic
    if h < g/l
        theta_max = acos(-h*l/g);
        k = sin(theta_max/2); % Solution as if maximum angle were initial angle
        e0 = ellipticK(k^2) - ellipticF(asin(sin(alpha0/2)/k), k^2); 
        [~,cn,dn] = ellipj(-sgn*sqrt(g/l)*t+e0, k^2); % Expression evaluated at (t-t0)
        cd = cn./dn;
        alpha = 2*asin(k*cd);
    % Bounded
    elseif h == g/l 
        e0 = atanh(sin(alpha0/2));
        alpha = 2*asin(tanh(sgn*sqrt(g/l)*t+e0)); % Expression evaluated at (t-t0)
    % Unbounded
    else
        m = 2*g/l/(h+g/l);
        e0 = ellipticF(sin(alpha0/2),m);
        [sn,cn,~] = ellipj(sgn*sqrt((h+g/l)/2)*t+e0, m); % Expression evaluated at (t-t0)
        alpha = 2*sign(cn).*asin(sn);
    end
end