function [val] = initial_theta_0(error, M, N, dt, theta)
    % 13/03/2020 18:10
    val = 0;
    for i = 1:M
        for j = 1:N
            val = val + (error(i,j+1)-(error(i,j)-theta*dt))^2;
        end
    end

end
