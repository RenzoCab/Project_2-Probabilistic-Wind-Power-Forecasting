function [Table_Training, new_bat] = new_batch_fixed(Table_Training, batch_size, N)
    % 11/02/2020 09:35

    Forecast     = Table_Training.Forecast;
    Forecast_Dot = Table_Training.Forecast_Dot;
    Error        = Table_Training.Error;
    new_bat      = [];
        
    for i = 1:batch_size
        for j = 2:N
            forecast(2*j-3:2*j-2)     = Forecast(i,j-1:j);
            forecast_dot(2*j-3:2*j-2) = Forecast_Dot(i,j-1:j);
            error(2*j-3:2*j-2)        = Error(i,j-1:j);
        end
        new_bat = [new_bat, [forecast; forecast_dot; error]];
    end
end