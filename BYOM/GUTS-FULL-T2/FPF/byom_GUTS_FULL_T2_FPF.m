%% BYOM, byom_guts_timevar_diazinon.m, full GUTS model
%
% * Author: Tjalling Jager 
% * Date: September 2020
% * Web support: <http://www.debtox.info/byom.html>
% * Back to index <walkthrough_guts.html>
% 
% MDA Changes: 
% Date: Mai 2021
% Model extended to use datasets at different temperatures (see derivatives).
%
% BYOM is a General framework for simulating model systems. The files in
% this directory use the ODE solver only, and therefore <simplefun.html
% simplefun.m> will be missing.
%
% *The model:* fitting survival data with the
% <http://www.debtox.info/about_guts.html GUTS> special cases based on the
% full model (TK and damage dynamics separate): SD, IT and mixed (or GUTS
% proper). 
%
%  Copyright (c) 2012-2020, Tjalling Jager, all rights reserved.
%  This source code is licensed under the MIT-style license found in the
%  LICENSE.txt file in the root directory of BYOM. 


%% Initial things
% Make sure that this script is in a directory somewhere *below* the BYOM
% folder.

clear, clear global % clear the workspace and globals
global DATA W X0mat % make the data set and initial states global variables
global glo          % allow for global parameters in structure glo
global pri zvd      % global structures for optional priors and zero-variate data
diary off           % turn of the diary function (if it is accidentaly on)
set(0,'DefaultFigureWindowStyle','docked'); % collect all figure into one window with tab controls
% set(0,'DefaultFigureWindowStyle','normal'); % separate figure windows

pathdefine % set path to the BYOM/engine directory
glo.basenm  = mfilename; % remember the filename for THIS file for the plots
glo.saveplt = 0; % save all plots as (1) Matlab figures, (2) JPEG file or (3) PDF (see all_options.txt)

%% The data set
% Data are entered in matrix form, time in rows, scenarios (exposure
% concentrations) in columns. 

% * -1 for multinomial likelihood (for survival data)
% * 0  for log-transform the data, then normal likelihood
% * 0.5 for square-root transform the data, then normal likelihood
% * 1  for no transformation of the data, then normal likelihood

% observed number of survivors, time in days, different scenario concentrations 
% are defined in Cw1-Cw3
DATA{1,1} = [-1	0	1	2	3	4	5       
            0	55	33	33	33	33	33      
            2	55	33	33	33	32	33      
            7	55	33	33	32	31	33      
            9	55	33	33	32	31	32      
            14	54	33	33	32	31	23      
            16	52	33	33	32	31	19      
            18	52	33	32	32	31	16      
            21	49	33	32	32	31	16      
            23	49	33	32	32	30	16      
            25	49	33	32	32	29	14      
            28	49	33	32	32	29	13   ];

DATA{2,1} = [-1	6	7	8	9	10	11      
            0	55	33	33	33	33	33     
            2	55	33	33	33	33	33      
            7	55	32	33	33	33	31      
            9	55	31	33	33	33	27     
            14	54	30	32	32	30	20      
            16	54	30	32	32	27	16      
            18	54	30	31	32	27	15      
            21	54	29	31	32	25	15      
            23	54	29	31	32	25	15      
            25	54	29	30	32	24	12      
            28	54	29	30	32	22	10   ];

DATA{3,1} = [-1	12	13	14	15	16	17
            0	55	33	33	33	33	33
            2	55	33	33	33	33	32
            7	55	33	33	33	33	28
            9	55	33	32	32	33	25
            14	54	32	30	28	32	18
            16	51	32	30	28	31	12
            18	50	32	29	28	30	10
            21	50	31	29	26	28	10
            23	50	31	29	26	26	 9
            25	50	30	29	25	24	 8
            28	47	30	28	24	22	 6   ];

% Exposure scenarios (to be replaced with measured data)
% 7 degrees experiment
Cw1 = [1	0	1       2       3       4       5
       0	0	0.28	0.93	2.82	9.00	29.58
       7	0	0.32	0.97	3.17	10.24	33.12
      14	0	0.30	1.00	3.17	10.33	31.33
      21	0	0.27	0.92	2.92	9.64	29.88
      28	0	0.30	1.03	3.17	10.59	32.01 ];

% 11 degrees expeiment
Cw2 = [1	6	7       8       9       10      11
       0	0	0.32	1.04	3.10	10.36	32.02
       7	0	0.31	1.07	3.27	10.79	33.05
      14	0	0.30	1.11	3.12	10.58	31.81
      21	0	0.29	1.03	3.24	10.66	32.75
      28	0	0.31	1.01	2.98	9.98	32.31 ];

%15 degrees experiment
Cw3 = [1	12	13      14      15      16      17
       0	0	0.29	1.04	3.13	10.70	31.37
       7	0	0.30	1.03	3.19	10.82	33.43
      14	0	0.30	0.91	2.85	9.86	30.20
      21	0	0.31	1.03	3.18	10.65	32.78
      28	0	0.26	0.92	2.90	9.58	31.74 ];

make_scen(1,Cw1,Cw2,Cw3); 

% no internal concentrations measured, TK parameters used from other experiment
% DATA{1,2} = NA;
% DATA{2,2} = NA;
% DATA{3,2} = NA;

% scaled damage; there are no data for this stage, but dummy data are
% added to make sure that the damage prediction is plotted (without having
% to plot ALL model curves in all sub-plots)
DATA{1,3} = [1  0     1     2     3     4     5    
             0  1e-6  1e-6  1e-6  1e-6  1e-6  1e-6  ]; 

DATA{2,3} = [1  6     7     8     9     10    11    
             0  1e-6  1e-6  1e-6  1e-6  1e-6  1e-6  ]; 

DATA{3,3} = [1  12    13    14     15   16    17
             0  1e-6  1e-6  1e-6  1e-6  1e-6  1e-6  ]; 

glo.wts = [1 1 0
           1 1 0
           1 1 0];  % set zero weight to fake data for scaled damage

% Experimental temperatures for each scenario
glo.Temp_scen = [  0      1      2      3      4      5      6      7      8      9     10     11     12     13     14     15     16     17    
                 280.15 280.15 280.15 280.15 280.15 280.15 284.15 284.15 284.15 284.15 284.15 284.15 288.15 288.15 288.15 288.15 288.15 288.15 ];   


% Create a table with nicer labels for the legends
Scenario = [0;1;2;3;4;5;
            6;7;8;9;10;11;
            12;13;14;15;16;17];
        
Label = {'control at 7C';'0.3 ug/L at 7C';'1 ug/L at 7C';'3 ug/L at 7C';'10 ug/L at 7C';'30 ug/L at 7C';
         'control at 11C';'0.3 ug/L at 11C';'1 ug/L at 11C';'3 ug/L at 11C';'10 ug/L at 11C';'30 ug/L at 11C'; 
         'control at 15C';'0.3 ug/L at 15C';'1 ug/L at 15C';'3 ug/L at 15C';'10 ug/L at 15C';'30 ug/L at 15C'};
   
glo.LabelTable = table(Scenario,Label); % create a Matlab table for the labels

%% Initial values for the state variables
% Initial states, scenarios in columns, states in rows. First row are the
% 'names' of all scenarios.

X0mat(1,:) = [Scenario]; % scenarios (concentrations or identifiers)
X0mat(2,:) = 1;           % initial survival probability
X0mat(3,:) = 0;           % initial internal concentration
X0mat(4,:) = 0;           % initial scaled damage

%% Initial values for the model parameters
% Model parameters are part of a 'structure' for easy reference. 

% global parameters for GUTS purposes
glo.sel  = 1; % select death mechanism: 1) SD 2) IT 3) mixed
glo.locS = 1 ; % location of survival probability in the state variable list
glo.locC = 2; % location of internal concentration in the state variable list
glo.locD = 3; % location of scaled damage in the state variable list
glo.fastrep = 0; % set to 1 to assume fast damage repair (death is driven by Ci)

% start values for SD
% syntax: par.name = [startvalue fit(0/1) minval maxval optional:log/normal scale (0/1)];
par.ku     = [     1.444  0      0.001         10 1]; 
par.ke     = [    0.1949  0      0.001         10 1]; 
par.T_A_tk = [      6102  0       1000      20000 1]; 
par.kr     = [     1.986  1      0.001         10 0]; 
par.mi     = [     143.8  1          0      1e+06 1]; 
par.hb     = [   0.01111  1          0      1e+06 1]; 
par.bi     = [  0.000692  1      1e-06      1e+06 0]; 
par.Fs     = [     38.55  1          1        100 1]; 
par.T_A_td = [     12010  1       1000      20000 1]; 

switch glo.sel % make sure that right parameters are fitted
    case 1 % SD
        par.Fs(2) = 0; % never fit the NEC spread!
    case 2 % IT
        par.bi(2) = 0; % never fit the killing rate!
    case 3 % mixed
        % do nothing: fit all parameters
end
if glo.fastrep == 1 % when damage repair is fast ...
    par.kr(2) = 0;  % never fit the repair rate!
end

%% Time vector and labels for plots
% Specify what to plot. If time vector glo.t is not specified, a default is
% constructed, based on the data set.

% specify the y-axis labels for each state variable
glo.ylab{1} = 'survival probability';
glo.ylab{2} = 'internal concentration (ug/kg)';
glo.ylab{3} = 'scaled damage (ug/kg)';
% specify the x-axis label (same for all states)
glo.xlab    = 'time (days)';
glo.leglab1 = 'scenario '; % legend label before the 'scenario' number
glo.leglab2 = 'ug/L'; % legend label after the 'scenario' number

prelim_checks % script to perform some preliminary checks and set things up
% Note: prelim_checks also fills all the options (opt_...) with defauls, so
% modify options after this call, if needed.

%% Calculations and plotting
% Here, the function is called that will do the calculation and the plotting.
% Options for the plotting can be set using opt_plot (see prelim_checks.m).
% Options for the optimsation routine can be set using opt_optim. Options
% for the ODE solver are part of the global glo. 
% 
% NOTE: for this package, the options useode and eventson in glo will not
% be functional: the ODE solver is always used, and the events function as
% well.

glo.stiff    = [0 1]; % ODE solver 0) ode45 (standard), 1) ode113 (moderately stiff), 2) ode15s (stiff)
% Seems that standard solver performs good enough, with tight tolerances.
% Second argument is for normally tight (1), tighter (2), or very tight (3)
% tolerances.
glo.break_time = 1; % break time vector up for ODE solver (1) or don't (0)

opt_optim.fit  = 1; % fit the parameters (1), or don't (0)
opt_optim.it   = 0; % show iterations of the optimisation (1, default) or not (0)
opt_plot.sho   = 0; % set to 1 to show all scenarios, 0 to only plot model for scenarios with data
opt_plot.annot = 1; % extra subplot in multiplot for fits: 1) box with parameter estimates, 2) overall legend

% optimise and plot (fitted parameters in par_out)
par_out = calc_optim(par,opt_optim); % start the optimisation
[Xall,Xall2x,Xall2y] = calc_and_plot(par_out,opt_plot); % calculate model lines and plot them

opt_tktd.max_exp = 0; % set to 1 to maximise exposure/damage plots on exposure rather than damage
plot_tktd(par_out,opt_tktd,[]); % leaving opt_conf empty suppresses all CIs for these plots

%% Likelihood region
% Make intervals on model predictions by using a sample of parameter sets
% taken from the joint likelihood-based conf. region. 

% % % UNCOMMENT FOLLOWING LINE(S) TO CALCULATE
% opt_likreg.subopt = 10; % number of sub-optimisations to perform to increase robustness
% opt_likreg.detail = 2; % detailed (1) or a coarse (2) calculation
% opt_likreg.burst  = 1000; % number of random samples from parameter space taken every iteration
% opt_likreg.lim_out = 1; % set to 1 to sample from a smaller part of space (enough for forward predictions)
% % 
% calc_likregion(par_out,1000,opt_likreg); % second argument number of samples (-1 to re-use saved sample from previous runs)

%% Plot results with confidence intervals
% The following code can be used to make plots with confidence intervals.
% Options for confidence bounds on model curves can be set using opt_conf
% (see prelim_checks). The plot_tktd function makes multiplots for the
% survival data, which are more readable when plotting with various
% intervals.

% % % UNCOMMENT FOLLOWING LINE(S) TO CALCULATE AND PLOT
% opt_conf.type    = 2; % make intervals from 1) slice sampler, 2)likelihood region, 3) parspace explorer
% opt_conf.lim_set = 2; % for lik-region sample: use limited set of n_lim points (1) or outer hull (2) to create CIs
% opt_tktd.max_exp = 0; % set to 1 to maximise exposure/damage plots on exposure rather than damage
% % 
% plot_tktd(par_out,opt_tktd,opt_conf); 


% %% Calculate LCx versus time
% % Here, the LCx (by default the LC50) is calculated at several time points.
% % LCx values are also printed on screen. If a sample from parameter space
% % is available (e.g., from the slice sampler or the likelihood region), it
% % can be used to calculate confidence bounds. 
% % 
% % Options for LCx (with confidence bounds) can be set using opt_ecx (see
% % prelim_checks). Note that opt_conf.type=-1 skips CIs.
% 
% opt_conf.type    = 2; % make intervals from 1) slice sampler, 2)likelihood region, 3) parspace explorer
% opt_conf.lim_set = 2; % for lik-region sample: use limited set of n_lim points (1) or outer hull (2) to create CIs
% opt_lcx_lim.Feff = 0.50; % effect level (>0 en <1), x/100 in LCx
% 
% % This is the general method as the fast methods won't work for the full
% % model (or at least: won't be faster than the general method).
% opt_ecx.Feff      = [0.50]; % effect levels (>0 en <1), x/100 in ECx
% opt_ecx.notitle   = 1; % set to 1 to suppress titles above ECx plots
% Tend = [2:0.2:3 3.5 4:8 10 15 22 28 40]; % times at which to calculate LCx, relative to control
% 
% calc_ecx(par_out,Tend,opt_ecx,opt_conf); % general method for ECx values
