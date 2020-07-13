function [est] = initial_guess(real_prod, M, N, dt)
    % 09/02/2020 09:30
    est = 0;

    for i = 1:M
        numerator = 0; denominator = 0;
        for j = 1:N
            numerator   = numerator + (real_prod(i,j+1) - real_prod(i,j))^2;
            denominator = denominator + real_prod(i,j)*(1-real_prod(i,j));
        end
        est = est + numerator/denominator;
    end
    est = est / (2*M*dt);

end
