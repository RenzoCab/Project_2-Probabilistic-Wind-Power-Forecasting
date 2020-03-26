function [bad_days] = remove_days_eps(Forecast, eps)
    % 14/03/2020 16:20

    [M, N_ini] = size(Forecast);
    bad_days   = [];
        
    for i = 1:M
        counted = 0; % To add only 1 time the day.
        for j = 1:N_ini
            
            if counted == 0 && (Forecast(i,j) < eps || Forecast(i,j) > 1-eps)
                counted     = 1;
                bad_days(end+1) = i;
            end
            
        end
    end
end