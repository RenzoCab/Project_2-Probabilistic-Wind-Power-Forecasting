close all;
clear all;
clc;

[Ta_Tra_Comp, Ta_Test_Compl, Ta_Comp] = load_data();

Table_Training = Ta_Tra_Comp; % We copy it so we can modify it.

Date              = Ta_Tra_Comp.Date;
Time              = Ta_Tra_Comp.Time;
Forecast          = Ta_Tra_Comp.Forecast;
Forecast_Dot      = Ta_Tra_Comp.Forecast_Dot;
real_ADME         = Ta_Tra_Comp.Real_ADME;
Error             = Ta_Tra_Comp.Error;
Error_Transitions = Ta_Tra_Comp.Error_Transitions;
Lamparti_Data     = Ta_Tra_Comp.Error_Lamp;
Lamparti_Tran     = Ta_Tra_Comp.Error_Lamp_Transitions;

Forecast_vec = Forecast(:); 
Error_vec    = Error(:);

len            = 10;
mean_error     = zeros(len,1);
mean_abs_error = zeros(len,1);
counter        = zeros(len,1);
intervals      = [0:1/len:1];

for i = 1:length(Forecast_vec)
    
    for j = 1:len
    
        if intervals(j) <= Forecast_vec(i) &&  Forecast_vec(i) < intervals(j+1)
            mean_error(j)     = mean_error(j) + Error_vec(i);
            mean_abs_error(j) = mean_abs_error(j) + abs(Error_vec(i));
            counter(j)        = counter(j) + 1;
        end

    end
    
end

for i = 1:len
    
    mean_error(i)     = mean_error(i) / counter(i);
    mean_abs_error(i) = mean_abs_error(i) / counter(i);
    
end

bar([0.05:0.1:0.95],mean_error); grid minor;
title('Mean Error');
xlabel('Forecast value');
saveas(gcf,'../../Python/Represas_Data_2/Wind_Data/someResults/final/mean_error','epsc');

figure;
bar([0.05:0.1:0.95],mean_abs_error); grid minor;
title('Absolute Mean Error');
xlabel('Forecast value');
saveas(gcf,'../../Python/Represas_Data_2/Wind_Data/someResults/final/mean_abs_error','epsc');

figure('Renderer', 'painters', 'Position', [10 10 1200 600]);
plot(Forecast,Error,'.'); grid minor;
title('Error over Forecast');
xlabel('Forecast value'); ylabel('Error value');
saveas(gcf,'../../Python/Represas_Data_2/Wind_Data/someResults/final/error_over_forecast','epsc');