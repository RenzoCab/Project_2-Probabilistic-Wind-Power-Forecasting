function [val] = initial_theta(error, alpha, forecast, M, N, dt, theta)
    % 09/02/2020 09:58
    val = 0;
    for i = 1:M
        for j = 1:N
            p_dot   = (forecast(i,j+1) - forecast(i,j)) / dt;
            theta_T = theta_t(theta, alpha, forecast(i,j), p_dot);
            val     = val + (error(i,j+1)-error(i,j)+theta_T*error(i,j)*dt)^2;
        end
    end

end
