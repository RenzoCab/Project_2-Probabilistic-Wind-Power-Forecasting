function [minus_eps_samples] = create_samples_minus_eps(Forecast, Error, gamma)
    % 15/03/2020 21:07
    
    % This returns the error V for only the samples where the forecast is
    % inside [e,1-e].

    [M, N_ini]         = size(Forecast);
    minus_eps_samples  = -1*ones(M*3,N_ini);
    eps_forecast       = -1*ones(M*3,N_ini); % For testing.
    
    flag         = 0;
    counter1     = 0;
    max_count_2  = 0;
    only_one     = 0; 
    % (^^) This is to remove the row, in case it has only one value.
        
    for i = 1:M
        counter2 = 1;
        for j = 1:N_ini
            
            if (Forecast(i,j) >= gamma && Forecast(i,j) <= 1-gamma)
                if flag == 0
                    flag     = 1;
                    counter1 = counter1 + 1;
                    only_one = 1;
                else
                    only_one = 0;
                end
                minus_eps_samples(counter1,counter2) = Error(i,j);
                eps_forecast(counter1,counter2)      = Forecast(i,j);
                counter2                             = counter2 + 1;
                max_count_2                          = max(max_count_2,counter2);
            else
                flag     = 0;
                counter2 = 1;
                if only_one == 1
                    counter1 = counter1 - 1;
                    only_one = 0;
                end
            end
            
        end
    end
    
    minus_eps_samples  = minus_eps_samples(1:counter1, 1:max_count_2-1);
    eps_forecast       = eps_forecast(1:counter1, 1:max_count_2-1);
end