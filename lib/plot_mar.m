function plot_mar(fnev)

data=load(fnev);

figure
plot(data(:,1),data(:,3)), grid on