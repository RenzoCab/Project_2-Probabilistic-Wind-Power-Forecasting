function [val,accum] = value_function_QMM_eps(error, theta_t, dt) % 15/03/2020 19:13
    [M, N_ini] = size(error);
    val        = 0;
    accum      = 0;
    for i = 1:M
        for j = 1:N_ini-1
            if error(i,j+1) ~= -1
                val   = val + (error(i,j+1) - error(i,j)*exp(-theta_t*dt))^2;
                accum = accum + 1;
            end
        end
    end
    val = val/accum;
end