function [val, accum] = mean_regression_eps(error, dt) % 16/03/2020 11:08
    [M, N_ini]  = size(error);
    num   = 0; den = 0; 
    accum = 0;
    for i = 1:M
        for j = 1:N_ini-1
            if error(i,j+1) ~= -1
                num   = num + error(i,j) * (error(i,j) - error(i,j+1));
                den   = den + error(i,j)^2;
                accum = accum + 1;
            end
        end
    end
    val = num/(den*dt);
end
