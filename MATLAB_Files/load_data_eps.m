function [Table_Training_Complete] = load_data_eps(ep)
    
    % 18/03/2020 10:41

    [Table_Training_Complete, ~, ~] = load_data();

    Time      = Table_Training_Complete.Time;
    Forecast  = Table_Training_Complete.Forecast;
    Real_ADME = Table_Training_Complete.Real_ADME;

    dt         = Time(1,2);
    [M, ~] = size(Forecast);
    
    for i = 1:M
        Forecast(i,:)          = pe(Forecast(i,:), ep);
        Forecast_Dot(i,:)      = pe_dot(Forecast(i,:), dt);
        Error(i,:)             = Real_ADME(i,:) - Forecast(i,:);
        Error_Transitions(i,:) = Error(i,2:end) - Error(i,1:end-1);
    end

    Table_Training_Complete.Forecast          = Forecast;
    Table_Training_Complete.Forecast_Dot      = Forecast_Dot;
    Table_Training_Complete.Error             = Error;
    Table_Training_Complete.Error_Transitions = Error_Transitions;
    
end