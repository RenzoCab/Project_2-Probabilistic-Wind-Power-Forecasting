close all;
clear all;
clc;

v = [-1:0.01:1];
x = [0:1:200000];

for i = 1:length(v)
    for j = 1:length(x)
        sur(i,j) = log_dist(v(i),x(j),x(j));
    end
end

[X,Y] = meshgrid(v,x);
contourf(X,Y,sur',100); colorbar;

%% 1D Plot:

v = [-1:0.01:1];

for i = 1:length(v)
	plo(i) = dist(v(i),1,1);
end

plot(v,plo);

%% 1D Plot:

v = [-1:0.01:1];

for i = 1:length(v)
	plo(i) = log_dist(v(i),500,500);
end

plot(v,plo);