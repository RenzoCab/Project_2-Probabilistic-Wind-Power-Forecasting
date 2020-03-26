function [val] = mean_regression(error, M, N, dt)
    % 13/03/2020 18:30
    val = 0;
    for i = 1:M
        num = 0; den = 0;
        for j = 1:N
            num = num + error(i,j) * (error(i,j) - error(i,j+1));
            den = den + error(i,j)^2;
        end
        val = val + num/den;
    end
    val = val/(M*dt);
end
