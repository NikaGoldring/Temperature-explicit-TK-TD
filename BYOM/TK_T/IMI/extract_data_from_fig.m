clear all;
close all;
clc;

%% Specify file for data extraction
open('20220726_TK_T_IMI_24.fig') ; % choose corresponding file name of MATLAB figure
h = findobj(gcf, 'Type', 'line');

% NOTE: As there were no input data for the metabolite provided at 07
% degrees, the locations of the data for the plotted lines within the
% matlab figure are different. Make sure to use the correct locations by
% identifying them first (i.e., open xdata and ydata, first lines in
% each section)

%% Extract xdata
%xdata= get(h,'XData'); % open variables stored in XData for identification
% xdata1= get(h(1,1),'XData'); % for 07 degrees: line for XData (is the same for all lines, so only exported once) 
xdata1= get(h(4,1),'XData'); % for 18 and 24 degrees: line for XData (is the same for all lines, so only exported once) 

%% Extract ydata
%ydata= get(h,'YData'); % open variables stored in YData for identification
% Parent dataset
% %07 degrees
% ydata1= get(h(7,1),'YData'); % line of upper CI - parent
% ydata2= get(h(8,1),'YData'); % line of lower CI - parent
% ydata3= get(h(9,1),'YData'); % line of model fit - parent
%18 and 24 degrees
ydata1= get(h(10,1),'YData'); % line of upper CI - parent
ydata2= get(h(11,1),'YData'); % line of lower CI - parent
ydata3= get(h(12,1),'YData'); % line of model fit - parent

%Metabolite dataset
% %07 degrees
% ydata4= get(h(1,1),'YData'); % line of upper CI - metabolite
% ydata5= get(h(2,1),'YData'); % line of lower CI - metabolite
% ydata6= get(h(3,1),'YData'); % line of model fit - metabolite
%18 and 24 degrees 
ydata4= get(h(4,1),'YData'); % line of upper CI - metabolite
ydata5= get(h(5,1),'YData'); % line of lower CI - metabolite
ydata6= get(h(6,1),'YData'); % line of model fit - metabolite

%% Export data to txt
fig01 = []; %create empty table
fig01(:,1) = xdata1 ; % xdata for all lines
fig01(:,2) = ydata1 ; % line of upper CI - parent
fig01(:,3) = ydata2 ; % line of lower CI - parent
fig01(:,4) = ydata3 ; % line of model fit - parent
fig01(:,5) = ydata4 ; % line of upper CI - metabolite
fig01(:,6) = ydata5 ; % line of lower CI - metabolite
fig01(:,7) = ydata6 ; % line of model fit - metabolite

dlmwrite('20220726_TK_T_IMI_24.txt', fig01, ','); % write dataframe in txt
