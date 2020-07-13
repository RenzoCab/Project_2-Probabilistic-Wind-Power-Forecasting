function [Table_Testing_Complete] = load_data_eps_test(ep,dataSet,dataKind,set_type)
    
    % 18/03/2020 10:41

    if strcmp(dataKind,'classic')
        if strcmp(set_type,'testing')
            [~, Table_Testing_Complete, ~] = load_data(dataSet);
        elseif strcmp(set_type,'training')
            [Table_Testing_Complete, ~, ~] = load_data(dataSet);
        end
    elseif strcmp(dataKind,'comparable')
        [~, Table_Testing_Complete, ~] = load_data_comparable(dataSet);
    end

    Time      = Table_Testing_Complete.Time;
    Forecast  = Table_Testing_Complete.Forecast;
    Real_ADME = Table_Testing_Complete.Real_ADME;

    dt         = Time(1,2);
    [M, ~] = size(Forecast);
    
    for i = 1:M
        Forecast(i,:)          = pe(Forecast(i,:), ep);
        Forecast_Dot(i,:)      = pe_dot(Forecast(i,:), dt);
        Error(i,:)             = Real_ADME(i,:) - Forecast(i,:);
        Error_Transitions(i,:) = Error(i,2:end) - Error(i,1:end-1);
    end

    Table_Testing_Complete.Forecast          = Forecast;
    Table_Testing_Complete.Forecast_Dot      = Forecast_Dot;
    Table_Testing_Complete.Error             = Error;
    Table_Testing_Complete.Error_Transitions = Error_Transitions;
    
end