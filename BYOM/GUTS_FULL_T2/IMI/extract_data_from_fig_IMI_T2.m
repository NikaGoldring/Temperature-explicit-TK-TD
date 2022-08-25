clear all;
close all;

%% Specify file for data extraction
open('20220808_GUTS_T2_SD_IMI_07.fig') ; % choose corresponding file name of MATLAB figure
h = findobj(gcf, 'Type', 'line'); %Get model line data
q = findobj(gcf, 'Type', 'ErrorBar');%Get data points and error bars

%% Extract xdata
%xdata= get(h,'XData'); % open variables stored in XData for identification
xdata= get(h(1,1),'XData'); % time points for model lines
Dxdata= get(q(1,1),'XData'); % time points for measured data

%% Extract ydata
% Lines
%ydata= get(h,'YData');
ydata2= get(h(1,1),'YData'); % 30 line 
ydata3= get(h(2,1),'YData'); % 30 upper CI
ydata4= get(h(3,1),'YData'); % 30 lower CI

ydata5= get(h(10,1),'YData'); % 10 line
ydata6= get(h(11,1),'YData'); % 10 upper CI
ydata7= get(h(12,1),'YData'); % 10 lower CI

ydata8= get(h(19,1),'YData'); % 3 line 
ydata9= get(h(20,1),'YData'); % 3 upper CI
ydata10= get(h(21,1),'YData'); % 3 lower CI

ydata11= get(h(28,1),'YData'); % 1 line
ydata12= get(h(29,1),'YData'); % 1 upper CI
ydata13= get(h(30,1),'YData'); % 1 lower CI

ydata14= get(h(37,1),'YData'); % 3 line 
ydata15= get(h(38,1),'YData'); % 3 upper CI
ydata16= get(h(39,1),'YData'); % 3 lower CI

ydata17= get(h(46,1),'YData'); % c line 
ydata18= get(h(47,1),'YData'); % c upper CI
ydata19= get(h(48,1),'YData'); % c lower CI


% Data and error bars
Dydata30= get(q(1,1),'YData'); % 30 
Dydata10= get(q(2,1),'YData'); % 10
Dydata3= get(q(3,1),'YData'); % 3
Dydata1= get(q(4,1),'YData'); % 1
Dydata03= get(q(5,1),'YData'); % 0.3
DydataC= get(q(6,1),'YData'); % control
% Negative error 
DNydata30= get(q(1,1),'YNegativeDelta'); % 30 
DNydata10= get(q(2,1),'YNegativeDelta'); % 10
DNydata3= get(q(3,1),'YNegativeDelta'); % 3
DNydata1= get(q(4,1),'YNegativeDelta'); % 1
DNydata03= get(q(5,1),'YNegativeDelta'); % 0.3
DNydataC= get(q(6,1),'YNegativeDelta'); % control
% Positive error
DPydata30= get(q(1,1),'YPositiveDelta'); % 30 
DPydata10= get(q(2,1),'YPositiveDelta'); % 10
DPydata3= get(q(3,1),'YPositiveDelta'); % 3
DPydata1= get(q(4,1),'YPositiveDelta'); % 1
DPydata03= get(q(5,1),'YPositiveDelta'); % 0.3
DPydataC= get(q(6,1),'YPositiveDelta'); % control

%% Export data to txt
fig01 = []; %create empty table
fig01(:,1) = xdata ; % xdata for all lines
% lines
fig01(:,2) = ydata2 ; % line
fig01(:,3) = ydata3 ; % upper CI
fig01(:,4) = ydata4 ; % lower CI
fig01(:,5) = ydata5 ; % line
fig01(:,6) = ydata6 ; % upper CI
fig01(:,7) = ydata7 ; % lower CI
fig01(:,8) = ydata8 ; % line
fig01(:,9) = ydata9 ; % upper CI
fig01(:,10) = ydata10 ; % lower CI
fig01(:,11) = ydata11 ; % line
fig01(:,12) = ydata12 ; % upper CI
fig01(:,13) = ydata13 ; % lower CI
fig01(:,14) = ydata14 ; % line
fig01(:,15) = ydata15 ; % upper CI
fig01(:,16) = ydata16 ; % lower CI
fig01(:,17) = ydata17 ; % line
fig01(:,18) = ydata18 ; % upper CI
fig01(:,19) = ydata19 ; % lower CI

%data and error bars
fig02 = []; %create empty table
fig02(:,1) = Dxdata ; % xdata for all model lines

fig02(:,2) = Dydata30 ; 
fig02(:,3) = Dydata10 ;
fig02(:,4) = Dydata3 ;
fig02(:,5) = Dydata1 ;
fig02(:,6) = Dydata03 ;
fig02(:,7) = DydataC ; 
% Negative error 
fig02(:,8) = DNydata30; % 30 
fig02(:,9) = DNydata10; % 10
fig02(:,10) = DNydata3; % 3
fig02(:,11) = DNydata1; % 1
fig02(:,12) = DNydata03; % 0.3
fig02(:,13) = DNydataC; % control
% Positive error
fig02(:,14) = DPydata30; % 30 
fig02(:,15) = DPydata10; % 10
fig02(:,16) = DPydata3; % 3
fig02(:,17) = DPydata1; % 1
fig02(:,18) = DPydata03; % 0.3
fig02(:,19) = DPydataC; % control
% 
dlmwrite('20220822_IMI_GUTS_SD_T2_model_07.txt', fig01, ','); %write dataframe in txt
dlmwrite(['20220822_IMI_GUTS_SD_T2_data_07.txt'], fig02, ','); %write dataframe in txt