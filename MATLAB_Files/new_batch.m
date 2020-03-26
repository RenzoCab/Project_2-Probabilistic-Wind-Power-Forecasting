function [Table_Training, new_bat] = new_batch(Table_Training, batch_size, N)
    % 03/02/2020 17:41
    idx          = randperm(height(Table_Training), batch_size); % Sample indices.
    idx          = -sort(-idx); % We order the indices from large to small.
    Forecast     = Table_Training.Forecast;
    Forecast_Dot = Table_Training.Forecast_Dot;
    Error        = Table_Training.Error;
    new_bat      = [];
        
    for i = 1:length(idx)
        Table_Training(idx(i),:) = []; % We remove the row that we sample from.
        for j = 2:N
            forecast(2*j-3:2*j-2)     = Forecast(idx(i),j-1:j);
            forecast_dot(2*j-3:2*j-2) = Forecast_Dot(idx(i),j-1:j);
            error(2*j-3:2*j-2)        = Error(idx(i),j-1:j);
        end
        new_bat = [new_bat, [forecast; forecast_dot; error]];
    end
end