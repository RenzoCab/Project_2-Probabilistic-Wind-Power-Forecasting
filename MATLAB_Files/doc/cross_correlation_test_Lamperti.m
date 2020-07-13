close all;
clear all;
clc;

% dataSet can be AWSTP (B), MTLOG (A) or UTEP5 (C).
dataSet       = 'MTLOG';
% epsilon can be 0.035 (B), 0.020 (A) or 0.028 (C).
epsilon       = 0.020;
% dataKind can be classic or comparable.
dataKind      = 'comparable';
train_data             = load_data_eps(epsilon,dataSet,dataKind);
test_data              = load_data_eps_test(epsilon,dataSet,dataKind);
[~, ~, table_complete] = load_data(dataSet);
lag                    = 150;
font_size              = 14;

train_error = transpose(train_data.Error_Lamp);
test_error  = transpose(test_data.Error_Lamp);
all_error   = transpose(table_complete.Error_Lamp);

[N1,M1]   = size(train_error);
[N2,M2]   = size(test_error);
all_error = all_error(:,end-(M1+M2)+1:end);
[N3,M3]   = size(all_error);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:M1-1
    train_xcorr(:,i) = xcorr(train_error(:,i),train_error(:,i+1),lag,'normalized');
end
for i = 1:M2-1
    test_xcorr(:,i) = xcorr(test_error(:,i),test_error(:,i+1),lag,'normalized');
end
for i = 1:M3-1
    all_xcorr(:,i) = xcorr(all_error(:,i),all_error(:,i+1),lag,'normalized');
end

X1 = 1:M1-1; X2 = 1:M2-1; X3 = 1:M3-1;
Y  = -lag:1:lag;
[x1,y1] = meshgrid(X1,Y);
[x2,y2] = meshgrid(X2,Y);
[x3,y3] = meshgrid(X3,Y);

figure('Renderer', 'painters', 'Position', [10 10 900 600]);
contour(x1,y1,train_xcorr); colorbar;
figure('Renderer', 'painters', 'Position', [10 10 900 600]);
imagesc(X1,Y,train_xcorr); colorbar;
title('Sample xcorr Function (observatios - training)');
xlabel('Day'); ylabel('Lag');
set(gca,'FontSize',font_size);
saveas(gcf,[pwd '/Results/cross_correlation/obs_training_L'],'epsc');

figure('Renderer', 'painters', 'Position', [10 10 900 600]);
contour(x2,y2,test_xcorr); colorbar;
figure('Renderer', 'painters', 'Position', [10 10 900 600]);
imagesc(X2,Y,test_xcorr); colorbar;
title('Sample xcorr Function (observatios - testing)');
xlabel('Day'); ylabel('Lag');
set(gca,'FontSize',font_size);
saveas(gcf,[pwd '/Results/cross_correlation/obs_testing_L'],'epsc');

figure('Renderer', 'painters', 'Position', [10 10 900 600]);
contour(x3,y3,all_xcorr); colorbar;
figure('Renderer', 'painters', 'Position', [10 10 900 600]);
imagesc(X3,Y,all_xcorr); colorbar;
title('Sample xcorr Function (observatios - all data)');
xlabel('Day'); ylabel('Lag');
set(gca,'FontSize',font_size);
saveas(gcf,[pwd '/Results/cross_correlation/obs_alldata_L'],'epsc');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

train_tran = transpose(train_data.Error_Lamp_Transitions);
test_tran  = transpose(test_data.Error_Lamp_Transitions);
all_tran   = transpose(table_complete.Error_Lamp_Transitions);

[N1,M1]   = size(train_tran);
[N2,M2]   = size(test_tran);
all_error = all_tran(:,end-(M1+M2)+1:end);
[N3,M3]   = size(all_tran);

for i = 1:M1-1
    train_xcorr(:,i) = xcorr(train_tran(:,i),train_tran(:,i+1),lag,'normalized');
end
for i = 1:M2-1
    test_xcorr(:,i) = xcorr(test_tran(:,i),test_tran(:,i+1),lag,'normalized');
end
for i = 1:M3-1
    all_xcorr(:,i) = xcorr(all_tran(:,i),all_tran(:,i+1),lag,'normalized');
end

X1 = 1:M1-1; X2 = 1:M2-1; X3 = 1:M3-1;
Y  = -lag:1:lag;
[x1,y1] = meshgrid(X1,Y);
[x2,y2] = meshgrid(X2,Y);
[x3,y3] = meshgrid(X3,Y);

figure('Renderer', 'painters', 'Position', [10 10 900 600]);
contour(x1,y1,train_xcorr); colorbar;
figure('Renderer', 'painters', 'Position', [10 10 900 600]);
imagesc(X1,Y,train_xcorr); colorbar;
title('Sample xcorr Function (transitions - training)');
xlabel('Day'); ylabel('Lag');
set(gca,'FontSize',font_size);
saveas(gcf,[pwd '/Results/cross_correlation/tran_training_L'],'epsc');

figure('Renderer', 'painters', 'Position', [10 10 900 600]);
contour(x2,y2,test_xcorr); colorbar;
figure('Renderer', 'painters', 'Position', [10 10 900 600]);
imagesc(X2,Y,test_xcorr); colorbar;
title('Sample xcorr Function (transitions - testing)');
xlabel('Day'); ylabel('Lag');
set(gca,'FontSize',font_size);
saveas(gcf,[pwd '/Results/cross_correlation/tran_testing_L'],'epsc');

figure('Renderer', 'painters', 'Position', [10 10 900 600]);
contour(x3,y3,all_xcorr); colorbar;
figure('Renderer', 'painters', 'Position', [10 10 900 600]);
imagesc(X3,Y,all_xcorr); colorbar;
title('Sample xcorr Function (transitions - all data)');
xlabel('Day'); ylabel('Lag');
set(gca,'FontSize',font_size);
saveas(gcf,[pwd '/Results/cross_correlation/tran_alldata_L'],'epsc');