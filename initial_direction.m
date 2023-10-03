function sgn = initial_direction(d_alpha0,alpha0,alpha_min,alpha_max)
% Function determines sign of initial polar displacement
% Used for all formulas with +-
    % With initial velocity
    if d_alpha0
        sgn = sign(d_alpha0);
    % At left turning point
    elseif abs(alpha0-alpha_min) < abs(alpha0-alpha_max)
        sgn = 1;
    % At right turning point
    else
        sgn = -1;
    end
end