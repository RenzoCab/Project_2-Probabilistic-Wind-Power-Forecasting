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

train_error = transpose(train_data.Error);
test_error  = transpose(test_data.Error);
all_error   = transpose(table_complete.Error);

train_error_vec = train_error(:);
test_error_vec  = test_error(:);
all_error_vec   = all_error(:);
L               = length(train_error_vec) + length(test_error_vec);
all_error_vec   = all_error_vec(end-L+1:end);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure;

[acf,lags,bounds] = autocorr(train_error_vec,'NumLags',200);
stem(lags,acf); xlabel('Lag'); ylabel('\rho(k)');
hold on;
% h = line(lags,bounds(1)*ones(length(acf),1));
% h1 = line(lags,bounds(2)*ones(length(acf),1));
% set(h,'color',[1 0 0]);
% set(h1,'color',[1 0 0]);

hold on;

[acf,lags,bounds] = autocorr(test_error_vec,'NumLags',200);
stem(lags,acf); xlabel('Lag'); ylabel('\rho(k)');
% h = line(lags,bounds(1)*ones(length(acf),1));
% h1 = line(lags,bounds(2)*ones(length(acf),1));
% set(h,'color',[1 0 0]);
% set(h1,'color',[1 0 0]);

[acf,lags,bounds] = autocorr(all_error_vec,'NumLags',200);
stem(lags,acf); xlabel('Lag'); ylabel('\rho(k)');
% h = line(lags,bounds(1)*ones(length(acf),1));
% h1 = line(lags,bounds(2)*ones(length(acf),1));
% set(h,'color',[1 0 0]);
% set(h1,'color',[1 0 0]);

grid minor;
title('Sample Autocorrelation Function (observations)');
legend({'Training data','Testing data','All data'});
saveas(gcf,[pwd '/Results/Autocorrelation/observations'],'epsc');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

train_tran = transpose(train_data.Error_Transitions);
test_tran  = transpose(test_data.Error_Transitions);
all_tran   = transpose(table_complete.Error_Transitions);

train_tran_vec = train_tran(:);
test_tran_vec  = test_tran(:);
all_tran_vec   = all_tran(:);
all_tran_vec   = all_tran_vec(end-L+1:end);

figure;

[acf,lags,bounds] = autocorr(train_tran_vec,'NumLags',200);
stem(lags,acf); xlabel('Lag'); ylabel('\rho(k)');
hold on;
% h = line(lags,bounds(1)*ones(length(acf),1));
% h1 = line(lags,bounds(2)*ones(length(acf),1));
% set(h,'color',[1 0 0]);
% set(h1,'color',[1 0 0]);

hold on;

[acf,lags,bounds] = autocorr(test_tran_vec,'NumLags',200);
stem(lags,acf); xlabel('Lag'); ylabel('\rho(k)');
% h = line(lags,bounds(1)*ones(length(acf),1));
% h1 = line(lags,bounds(2)*ones(length(acf),1));
% set(h,'color',[1 0 0]);
% set(h1,'color',[1 0 0]);

[acf,lags,bounds] = autocorr(all_tran_vec,'NumLags',200);
stem(lags,acf); xlabel('Lag'); ylabel('\rho(k)');
% h = line(lags,bounds(1)*ones(length(acf),1));
% h1 = line(lags,bounds(2)*ones(length(acf),1));
% set(h,'color',[1 0 0]);
% set(h1,'color',[1 0 0]);

grid minor;
title('Sample Autocorrelation Function (transitions)');
legend({'Training data','Testing data','All data'});
saveas(gcf,[pwd '/Results/Autocorrelation/transitions'],'epsc');