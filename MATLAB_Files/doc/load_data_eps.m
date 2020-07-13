function [table_to_use] = load_data_eps(ep,dataSet,dataKind,dataNumber)
    
    % 19/04/2020 15:45
    
    if strcmp(dataKind,'classic')
        [Table_Training_Complete, Table_Testing_Complete, ~] = load_data(dataSet);
    elseif strcmp(dataKind,'comparable')
        [Table_Training_Complete, Table_Testing_Complete, ~] = load_data_comparable(dataSet);
    end
    
    if strcmp(dataNumber,'training')
        table_to_use = Table_Training_Complete;
    elseif strcmp(dataNumber,'testing')
        table_to_use = Table_Testing_Complete;
    end

    Time      = table_to_use.Time;
    Forecast  = table_to_use.Forecast;
    Real_ADME = table_to_use.Real_ADME;

    dt         = Time(1,2);
    [M, ~] = size(Forecast);
    
    for i = 1:M
        Forecast(i,:)          = pe(Forecast(i,:), ep);
        Forecast_Dot(i,:)      = pe_dot(Forecast(i,:), dt);
        Error(i,:)             = Real_ADME(i,:) - Forecast(i,:);
        Error_Transitions(i,:) = Error(i,2:end) - Error(i,1:end-1);
    end

    table_to_use.Forecast          = Forecast;
    table_to_use.Forecast_Dot      = Forecast_Dot;
    table_to_use.Error             = Error;
    table_to_use.Error_Transitions = Error_Transitions;
    
end