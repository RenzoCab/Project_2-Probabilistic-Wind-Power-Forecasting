close all;
clear all;
clc;

epsilon     = 0.018;
Ta_Tra_Comp = load_data_eps(epsilon);

Time              = Ta_Tra_Comp.Time;
Forecast_Dot      = Ta_Tra_Comp.Forecast_Dot;
Error             = Ta_Tra_Comp.Error;

Forecast_Transitions = Forecast_Dot(:) * Time(1,2);

h1 = histogram(Forecast_Transitions);
grid minor; title('Forecast Transitions');
h1.Normalization = 'pdf';
xlabel('Value of Transition');
ylabel('Probability');
saveas(gcf,[pwd '/Results/histograms/others/forecast_transitions'],'epsc');

figure;
h1 = histogram(Error);
grid minor; title('Error Measurement');
h1.Normalization = 'pdf';
xlabel('Value of Error');
ylabel('Probability');
saveas(gcf,[pwd '/Results/histograms/others/error_measurement'],'epsc');

Forecast_Transitions = Forecast_Transitions(:);
Error(:,end)         = []; % We remove the error at final time for each day.
Error                = Error(:);
figure;
plot(Forecast_Transitions,Error,'*');
xlim([-0.03 0.03]); ylim([-0.3 0.3]);
xlabel('Prev. Rel. Forecast'); ylabel('Prev. Error (V)');
grid minor;
saveas(gcf,[pwd '/Results/histograms/others/error_and_forecast'],'epsc');