load('MTLOG_0100_and_Real_24h.mat');

forecast  = Table.Forecast;
real_ute  = Table.Real_UTE;
real_adme = Table.Real_ADME;
time      = linspace(1,25,145);

figure(70);
currentDate = '20190101';
count = 1;

for i = 1:length(real_adme)
    
    set(0,'CurrentFigure',70); clf(70);
    P = plot(time,real_ute(i,:)); P.LineWidth = 1;
    hold on;
    P = plot(time,real_adme{i});  P.LineWidth = 1;
    P = plot(time,forecast(i,:)); P.LineWidth = 1;
    ylim([0 1400]); grid minor;
    xlabel('Hours (24 hrs starting at 01:00)'); ylabel('MW'); title(currentDate);
    xlim([1 25]);
    legend('Real UTE','Real ADME','Forecast');
    
%     saveas(gcf,[pwd '/someResults/dataToUse_24/',num2str(count)],'epsc');
    pause(0.1);
    
    currentDate = datetime(currentDate,'InputFormat','yyyyMMdd');
    currentDate = currentDate + days(1);
    currentDate = datestr(currentDate,'yyyymmdd');
    count = count + 1;
    
end