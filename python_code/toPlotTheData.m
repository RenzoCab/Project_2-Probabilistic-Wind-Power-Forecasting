close all
clear all
clc

data_p = readtable('data_p.xlsx');
data_X = readtable('data_X.xlsx');
data_V = readtable('data_V.xlsx');
div = 10;

% Ploting p:
figure;
hold on; title('data: p');
for i = 1:height(data_p)
    plot([1:59],table2array(data_p(i,:)));
end
figure;
hold on; title('data: p');
for i = 1:floor(height(data_p)/div)
    plot([1+(i-1)*59:59+(i-1)*59],table2array(data_p(i,:)));
end

% Ploting X:
figure;
hold on; title('data: X');
for i = 1:height(data_X)
    plot([1:59],table2array(data_X(i,:)));
end
figure;
hold on; title('data: X');
for i = 1:floor(height(data_X)/div)
    plot([1+(i-1)*59:59+(i-1)*59],table2array(data_X(i,:)));
end

% Ploting V:
figure;
hold on; title('data: V');
for i = 1:height(data_V)
    plot([1:59],table2array(data_V(i,:)));
end
figure;
hold on; title('data: V');
for i = 1:floor(height(data_V)/div)
    plot([1+(i-1)*59:59+(i-1)*59],table2array(data_V(i,:)));
end

%% Check if we have repeated data:

data_p = table2array(data_p);
data_X = table2array(data_X);
data_V = table2array(data_V);

unique_p = unique(data_p,'rows');
unique_X = unique(data_X,'rows');
unique_V = unique(data_V,'rows');