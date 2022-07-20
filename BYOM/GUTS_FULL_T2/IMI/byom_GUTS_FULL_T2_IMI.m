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
              2	54	33	33	33	33	32
              4	54	33	33	33	32	32
              7	54	33	33	33	32	32
              9	53	33	33	33	32	31
             11	53	33	33	33	32	30
             14	53	33	33	33	31	30
             16	52	33	33	33	31	28
             18	51	33	33	33	29	28
             21	51	33	32	33	29	28
             23	51	33	32	31	28	27
             28	49	33	32	29	25	26   ];

DATA{2,1} = [-1	6	7	8	9	10	11
              0	55	33	33	33	33	33
              2	55	33	33	33	33	30
              4	55	33	33	33	33	30
              7	55	33	33	33	32	30
              9	55	33	33	33	31	30
             11	55	33	33	33	30	29
             14	55	33	33	32	28	27
             16	55	33	32	32	27	26
             18	55	33	32	32	24	25
             21	55	33	31	32	22	22
             23	55	33	29	32	20	21
             28	54	33	29	30	15	19   ];

DATA{3,1} = [-1	12	13	14	15	16	17
              0	55	33	33	33	33	33
              2	55	32	33	33	33	31
              4	55	32	33	33	33	29
              7	54	31	32	31	32	29
              9	53	28	32	31	31	27
             11	53	27	32	31	31	26
             14	52	26	32	31	26	25
             16	52	26	32	28	26	23
             18	52	26	32	25	25	21
             21	51	24	31	25	23	17
             23	51	23	31	24	23	16
             28	50	23	25	21	20	13   ];

% Exposure scenarios (to be replaced with measured data)
% 7 degrees experiment
Cw1 = [1	0	1       2       3       4       5
       0	0	0.21	0.80	2.69	8.94	25.47
       7	0	0.20	0.73	2.35	8.61	23.44
      14	0	0.20	0.78	2.54	8.81	25.03
      21	0	0.30	0.99	2.89	9.27	27.24
      28	0	0.29	0.97	2.83	9.45	25.09  ];

% 11 degrees expeiment
Cw2 = [1	6	7       8       9       10      11
       0	0	0.22	0.85	2.76	9.21	24.18
       7	0	0.20	0.82	2.53	8.97	25.06
      14	0	0.21	0.84	2.57	9.07	25.42
      21	0	0.30	0.99	2.81	9.47	26.73
      28	0	0.31	1.01	2.84	9.47	26.84  ];

%15 degrees experiment
Cw3 = [1	12	13      14      15      16      17
       0	0	0.23	0.82	2.65	9.31	25.32
       7	0	0.22	0.80	2.49	9.13	24.99
      14	0	0.26	0.83	2.54	9.72	26.22
      21	0	0.29	0.98	2.81	9.30	25.51
      28	0	0.32	1.07	2.92	9.20	24.74  ];

make_scen(1,Cw1,Cw2,Cw3); 

% no internal concentrations measured, TK parameters used from other experiment
% DATA{1,2} = NA;
% DATA{2,2} = NA;
% DATA{3,2} = NA;

% scaled damage; there are no data for this stage, but dummy data are
% added to make sure that the damage prediction is plotted (without having
% to plot ALL model curves in all sub-plots)
DATA{1,3} = [1   0     1     2     3     4     5    
             0  1e-6  1e-6  1e-6  1e-6  1e-6  1e-6  ]; 

DATA{2,3} = [1   6     7     8     9     10    11    
             0  1e-6  1e-6  1e-6  1e-6  1e-6  1e-6  ]; 

DATA{3,3} = [1   12    13    14     15   16    17
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
glo.sel  = 2; % select death mechanism: 1) SD 2) IT 3) mixed
glo.locS = 1 ; % location of survival probability in the state variable list
glo.locC = 2; % location of internal concentration in the state variable list
glo.locD = 3; % location of scaled damage in the state variable list
glo.fastrep = 0; % set to 1 to assume fast damage repair (death is driven by Ci)

% start values for SD
% syntax: par.name = [startvalue fit(0/1) minval maxval optional:log/normal scale (0/1)];
par.ku    = [3.213    0 1e-3  10 1];  % uptake rate constant, L/kg/d
par.ke    = [0.1524  0 1e-3  10 1]; % elimination rate constant (d-1)
par.T_A_tk = [2818 0 1000 20000 1];  % Arrhenius temperature for TK rates, Kelvin (startvalue from AmP = 10000) 
par.kr     = [     0.001  1      0.001         10 0]; 
par.mi     = [    12.72  1          0      1e+06 1]; 
par.hb     = [ 0.01217 1          0      1e+06 1]; 
par.bi     = [     0.002  1      1e-06      1e+06 0]; 
par.Fs     = [       100.0  1          1        100 1]; 
par.T_A_td = [     1.178e+04  1       1000      20000 1]; 

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

% UNCOMMENT FOLLOWING LINE(S) TO CALCULATE
opt_likreg.subopt = 10; % number of sub-optimisations to perform to increase robustness
opt_likreg.detail = 2; % detailed (1) or a coarse (2) calculation
opt_likreg.burst  = 1000; % number of random samples from parameter space taken every iteration
opt_likreg.lim_out = 1; % set to 1 to sample from a smaller part of space (enough for forward predictions)
% 
calc_likregion(par_out,1000,opt_likreg); % second argument number of samples (-1 to re-use saved sample from previous runs)

%% Plot results with confidence intervals
% The following code can be used to make plots with confidence intervals.
% Options for confidence bounds on model curves can be set using opt_conf
% (see prelim_checks). The plot_tktd function makes multiplots for the
% survival data, which are more readable when plotting with various
% intervals.

% UNCOMMENT FOLLOWING LINE(S) TO CALCULATE AND PLOT
opt_conf.type    = 2; % make intervals from 1) slice sampler, 2)likelihood region, 3) parspace explorer
opt_conf.lim_set = 2; % for lik-region sample: use limited set of n_lim points (1) or outer hull (2) to create CIs
opt_tktd.max_exp = 0; % set to 1 to maximise exposure/damage plots on exposure rather than damage
% 
plot_tktd(par_out,opt_tktd,opt_conf); 
