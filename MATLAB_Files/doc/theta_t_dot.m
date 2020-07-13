function [theta_t_dot] = theta_t_dot(theta_0, p, p_dot, p_ddot)
    
    % 23/02/2020 01:01
    if theta_0 < abs(p_dot)/p && 2*p < 1
        theta_t_dot = (p*sign(p_dot)*p_ddot^2-abs(p_dot)*p_dot) / p^2;
    elseif theta_0 + abs(p_dot)/(p-1) < 0 && 2*p >= 1
        theta_t_dot = ((1-p)*sign(p_dot)*p_ddot^2-abs(p_dot)*p_dot) / (p-1)^2;
    else
        theta_t_dot = 0;
    end

end