close all
clear all
clc

time   = [2000:1:2019];
wind   = [17 24 31 39 48 59 74 94 121 159 196 239 283 319 370 423 489 537 600 666];
solar  = [1.3 1.6 2.1 2.6 3.7 5.1 6.7 9.2 16 23 40 70 101 138 178 290 285 340 400 467];
toPlot = [wind;solar];

figure('Renderer', 'painters', 'Position', [10 10 900 400]);
bar(time,toPlot);
grid minor;
legend('Wind','Solar','location','northwest');
ylabel('GW');
xlabel('Year');
title('Installed Wind and Solar Power over the years');

saveas(gcf,[pwd '/Results/wind_and_solar/plot_over_years'],'epsc');