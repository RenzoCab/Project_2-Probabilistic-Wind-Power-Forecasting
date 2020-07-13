function [Table_Training_Complete, Table_Testing_Complete, Table_Complete] = load_data(dataSet)
    
    % 19/04/2020 15:45
    % dataSet can be AWSTP, MTLOG or UTEP5.
    load(['../../Python/Represas_Data_2/Wind_Data/',dataSet,'_0100_and_Real_24h_Training_Data.mat']);
    load(['../../Python/Represas_Data_2/Wind_Data/',dataSet,'_0100_and_Real_24h_Testing_Data.mat']);
    load(['../../Python/Represas_Data_2/Wind_Data/',dataSet,'_0100_and_Real_24h_Complete_Data.mat']);

%     Date              = Table_Training_Complete.Date;
%     Time              = Table_Training_Complete.Time;
%     Forecast          = Table_Training_Complete.Forecast;
%     Forecast_Dot      = Table_Training_Complete.Forecast_Dot;
%     real_ADME         = Table_Training_Complete.Real_ADME;
%     Error             = Table_Training_Complete.Error;
%     Error_Transitions = Table_Training_Complete.Error_Transitions;
%     Lamparti_Data     = Table_Training_Complete.Error_Lamp;
%     Lamparti_Tran     = Table_Training_Complete.Error_Lamp_Transitions;

end