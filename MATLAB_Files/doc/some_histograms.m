close all;
clear all;
clc;

% dataSet can be AWSTP (B), MTLOG (A) or UTEP5 (C).
dataSet                = 'MTLOG';
% epsilon can be 0.035 (B), 0.018 (A) or 0.028 (C).
epsilon                = 0.018;
Ta_Tra_Comp            = load_data_eps(epsilon,dataSet);

Time              = Ta_Tra_Comp.Time;
Forecast_Dot      = Ta_Tra_Comp.Forecast_Dot;
Error             = Ta_Tra_Comp.Error;
Real_ADME         = Ta_Tra_Comp.Real_ADME;

Forecast_Transitions = Forecast_Dot(:) * Time(1,2);
Error_aux = Error;

% We create real measurements transitions:
[M,N] = size(Real_ADME);
for i = 1:M
    
    for j = 2:N
        
        Real_ADME_transitions(i,j-1) = Real_ADME(i,j) - Real_ADME(i,j-1);
        
    end
    
end

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

figure;
h1 = histogram(Real_ADME_transitions);
grid minor; title('Real Measurements Transitions');
h1.Normalization = 'pdf';
xlabel('Value of Transition');
ylabel('Probability');
saveas(gcf,[pwd '/Results/histograms/others/real_measurements_transitions'],'epsc');

figure;
h1 = histogram(Real_ADME);
grid minor; title('Real Measurements');
h1.Normalization = 'pdf';
xlabel('Value of Measurements');
ylabel('Probability');
saveas(gcf,[pwd '/Results/histograms/others/real_measurements'],'epsc');

figure;
h1 = histogram(Error_aux(:,1));
grid minor; title('Error Measurement Initial');
h1.Normalization = 'pdf';
xlabel('Value of Error');
ylabel('Probability');
saveas(gcf,[pwd '/Results/histograms/others/error_measurement_initial'],'epsc');