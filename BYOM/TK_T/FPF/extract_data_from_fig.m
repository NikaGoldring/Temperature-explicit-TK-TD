clear all;
close all;
clc;

%% Specify file for data extraction
open('20220726_TK_T_FPF_24.fig') ; % choose corresponding file name of MATLAB figure
h = findobj(gcf, 'Type', 'line');

%% Extract xdata
% xdata= get(h,'XData'); % open variables stored in XData for identification
xdata1= get(h(4,1),'XData'); % line for XData (is the same for all lines, so only exported once) 

%% Extract ydata
% ydata= get(h,'YData'); % open variables stored in YData for identification
ydata1= get(h(4,1),'YData'); % line of upper CI - parent
ydata2= get(h(5,1),'YData'); % line of lower CI - parent
ydata3= get(h(6,1),'YData'); % line of model fit - parent


%% Export data to txt
fig01 = []; %create empty table
fig01(:,1) = xdata1 ; % xdata for all lines
fig01(:,2) = ydata1 ; % line of upper CI - parent
fig01(:,3) = ydata2 ; % line of lower CI - parent
fig01(:,4) = ydata3 ; % line of model fit - parent

dlmwrite('20220726_TK_T_FPF_24.txt', fig01, ','); % write dataframe in txt
