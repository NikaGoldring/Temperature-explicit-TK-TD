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
% Date: 28.07.2022
% Model changed to be used by BYOMv.6.2 and added biotransformation 
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
% TK parameters (fixed)
par.ke     = [0.1397   0 1e-6 10 0];  % elimination rate constant, d-1
par.ku     = [3.250    0 0.01 10 0];  % uptake rate constant, L/kg/d
par.km     = [0.02172  0 1e-4 10 0];  % formation rate of the metabolite
%par.kem    = [0.4884   0 1e-4 10 0];  % elimination rate of the metabolite
par.T_A_tk = [3044     0 0.01 2e4 1];  % Arrhenius temperature, Kelvin
% % TD parameter starts for SD
% par.kr    = [0.069  1  1e-3  10    0]; % damage repair rate constant (d-1)
% par.mi    = [0.001  1  1e-3  20    0]; % median threshold for survival (nmol/kg)
% par.hb    = [0.004  1  1e-3  10    0]; % background hazard rate (d-1)
% par.bi    = [1e-3   1  1e-6  2e-3  0]; % killing rate (kg/nmol/d) (SD and mixed)
% par.Fs    = [42.99  0  1     100   1]; % fraction spread of threshold distribution (-) (IT and mixed)
% TD parameter starts for IT
par.kr    = [0.001    1 1e-3 1   0]; % damage repair rate constant (d-1)
par.mi    = [12.02    1 1    1e3 0]; % median threshold for survival (nmol/kg)
par.hb    = [0.003234 1 1e-6 1   0]; % background hazard rate (d-1)
par.bi    = [2e-3     1 1e-6 1e6 0]; % killing rate (kg/nmol/d) (SD and mixed)
par.Fs    = [42.19    1    1 1e3 1]; % fraction spread of threshold distribution (-) (IT and mixed)

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
% Here, the functions are called that will do the calculation and the
% plotting. 
% 
% Options for the plotting can be set using opt_plot (see prelim_checks.m).
% Options for the optimsation routine can be set using opt_optim. Options
% for the ODE solver are part of the global glo. For the demo, the
% iterations were turned off (opt_optim.it = 0).

opt_optim.type = 4; % optimisation method 1) simplex, 4) parameter-space explorer
opt_optim.fit  = 1; % fit the parameters (1), or don't (0)
opt_optim.it   = 0; % show iterations of the optimisation (1, default) or not (0)
glo.useode     = 1; % calculate model using ODE solver (1) or analytical solution (0)

opt_optim.ps_plots = 0; % when set to 1, makes intermediate plots to monitor progress of parameter-space explorer
opt_optim.ps_profs = 1; % when set to 1, makes profiles and additional sampling for parameter-space explorer
opt_optim.ps_rough = 1; % set to 1 for rough settings of parameter-space explorer, 0 for settings as in openGUTS
opt_optim.ps_saved = 0; % use saved set for parameter-space explorer (1) or not (0);

% optimise and plot (fitted parameters in par_out)
par_out = calc_optim(par,opt_optim); % start the optimisation
% no plotting here; we'll immediately plot with CIs below

%% Plot results with confidence intervals
% The following code can be used to make a standard plot (the same as for
% the fits), but with confidence intervals. Options for confidence bounds
% on model curves can be set using opt_conf (see prelim_checks).
% 
% Use opt_conf.type to tell calc_conf which sample to use: 
% -1) Skips CIs (zero does the same, and an empty opt_conf as well).
% 1) Bayesian MCMC sample (default); CI on predictions are 95% ranges on 
% the model curves from the sample 
% 2) parameter sets from a joint likelihood region using the shooting 
% method (limited sets can be used), which will yield (asymptotically) 95% 
% CIs on predictions
% 3) as option 2, but using the parameter-space explorer

opt_conf.type    = 3; % make intervals from 1) slice sampler, 2) likelihood region shooting, 3) parspace explorer
opt_conf.lim_set = 0; % use limited set of n_lim points (1) or outer hull (2, likelihood methods only) to create CIs
opt_conf.sens    = 0; % type of analysis 0) no sensitivities 1) corr. with state, 2) corr. with state/control, 3) corr. with relative change of state over time

out_conf = calc_conf(par_out,opt_conf);   % calculate confidence intervals on model curves
calc_and_plot(par_out,opt_plot,out_conf); % call the plotting routine again to plot fits with CIs

%% TKTD plots
opt_tktd.repls = 0;     % plot individual replicates (1) or means (0)
opt_tktd.min   = 0;     % set to 1 to show a dotted line for the control (lowest) treatment

plot_tktd(par_out,opt_tktd,opt_conf);

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
% 
% opt_ecx.id_sel    = [0   0 1]; % scenario to use from X0mat, scenario identifier, flag for ECx to use scenarios rather than concentrations
% % opt_ecx.id_sel    = [6   6 1]; % scenario to use from X0mat, scenario identifier, flag for ECx to use scenarios rather than concentrations
% % opt_ecx.id_sel    = [12 12 1]; % scenario to use from X0mat, scenario identifier, flag for ECx to use scenarios rather than concentrations
% 
% % This is the general method as the fast methods won't work for the full
% % model (or at least: won't be faster than the general method).
% opt_ecx.Feff      = [0.50]; % effect levels (>0 en <1), x/100 in ECx
% opt_ecx.notitle   = 1; % set to 1 to suppress titles above ECx plots
% Tend = [1:28]; % times at which to calculate LCx, relative to control
% 
% calc_ecx([],Tend,opt_ecx,opt_conf); % general method for ECx values


%% Other CI options
% % % %% Calculations and plotting
% % % % Here, the functions are called that will do the calculation and the
% % % % plotting. Note that calc_plot can provide all of the plotting information
% % % % as output, so you can also make your own customised plots. This section,
% % % % by default, makes a multiplot with all state variables (each in its own
% % % % panel of the multiplot). When the model is fitted to data, output is
% % % % provided on the screen (and added to the log-file results.out). The
% % % % zero-variate data point is also plotted with its prediction (although the
% % % % legend obscures it here).
% % % %
% % % % Options for the plotting can be set using opt_plot (see prelim_checks.m).
% % % % Options for the optimisation routine can be set using opt_optim. Options
% % % % for the ODE solver are part of the global glo.
% % % %
% % % % You can turn on the events function there too, to smoothly catch the
% % % % discontinuity in the model. For the demo, the iterations were turned off
% % % % (opt_optim.it = 0).
% % % 
% % % glo.eventson   = 1; % events function on (1) or off (0)
% % % glo.useode     = 1; % calculate model using ODE solver (1) or analytical solution (0)
% % % opt_optim.it   = 0; % show iterations of the optimisation (1, default) or not (0)
% % % opt_plot.annot = 1; % extra subplot in multiplot for fits: 1) box with parameter estimates, 2) overall legend
% % % opt_plot.bw    = 0; % if set to 1, plots in black and white with different plot symbols
% % % 
% % % % glo.diary = 'test_results.out'; % use a different name for the diary
% % % 
% % % % optimise and plot (fitted parameters in par_out)
% % % par_out = calc_optim(par,opt_optim); % start the optimisation
% % % calc_and_plot(par_out,opt_plot);     % calculate model lines and plot them
% % % 
% % % % save_plot(gcf,'fit_example',[],3) % save active figure as PDF (see save_plot for more options)
% % % 
% % % % return % stop here, and run analyses below manually later
% % % 
% % % %% Local sensitivity analysis
% % % % Local sensitivity analysis of the model. All model parameters are
% % % % increased one-by-one by a small percentage. The sensitivity score is
% % % % by default scaled (dX/X p/dp) or alternatively absolute (dX p/dp).
% % % %
% % % % Options for the sensitivity can be set using opt_sense (see
% % % % prelim_checks.m).
% % %  
% % % % % UNCOMMENT LINE(S) TO CALCULATE
% % % % opt_sens.state = 2; % vector with state variables for the sensitivity analysis (0 does all, if sens>0)
% % % % calc_localsens(par_out,opt_sens)
% % % 
% % % %% Profiling the likelihood
% % % % By profiling you make robust confidence intervals for one or more of your
% % % % parameters. Use the names of the parameters as they occurs in your
% % % % parameter structure _par_ above. This can be a single string (e.g.,
% % % % 'kd'), a cell array of strings (e.g., {'kd','ke'}), or 'all' to profile
% % % % all fitted parameters. This example produces a profile for each parameter
% % % % and provides the 95% confidence interval (on screen and indicated by the
% % % % horizontal broken line in the plot).
% % % %
% % % % Notes: profiling for complex models is a slow process, so grab a coffee!
% % % % If the profile finds a better solution, it breaks the analysis (as long
% % % % as you keep the default opt_prof.brkprof=1) and displays the parameters
% % % % for the new optimum on screen (and in results.out). For the NON-parallel
% % % % version of calc_proflik, the new optimum is also immediately saves to the
% % % % log file profiles_newopt.out. You can break of the analysis by pressing
% % % % ctrl-c anytime, and use the values from the log file to restart (copy the
% % % % better values into your _par_ structure). For parallel processing, saving
% % % % would need more thought.
% % % %
% % % % Options for profiling can be set using opt_prof (see prelim_checks.m).
% % % 
% % % opt_prof.detail   = 1; % detailed (1) or a coarse (2) calculation
% % % opt_prof.subopt   = 0; % number of sub-optimisations to perform to increase robustness
% % % opt_prof.brkprof  = 2; % when a better optimum is located, stop (1) or automatically refit (2)
% % % 
% % % % % UNCOMMENT LINE(S) TO CALCULATE
% % % % par_better = calc_proflik(par_out,'all',opt_prof,opt_optim);  % calculate a profile
% % % % if ~isempty(par_better)                 % if the profiling found a better optimum ...
% % % %     print_par(par_better) % display on screen in formatted manner, so it can be copied into the code
% % % %     calc_and_plot(par_better,opt_plot); % calculate model lines and plot them
% % % % end
% % % 
% % % %% Likelihood region
% % % % Another way to make intervals on model predictions is to use a sample of
% % % % parameter sets taken from the joint likelihood-based conf. region. This
% % % % is done by the function calc_likregion.m. It first does profiling of all
% % % % fitted parameters to find the edges of the region. Then, Latin-Hypercube
% % % % shooting, keeping only those parameter combinations that are not rejected
% % % % at the 95% level in a lik.-rat. test. The inner rim will be used for CIs
% % % % on forward predictions.
% % % %
% % % % Options for the likelihood region can be set using opt_likreg (see
% % % % prelim_checks.m). For the profiling part, use the options in opt_prof.
% % % 
% % % opt_prof.detail  = 1; % detailed (1) or a coarse (2) calculation
% % % opt_prof.subopt  = 0; % number of sub-optimisations to perform to increase robustness
% % % opt_prof.re_fit  = 1; % set to 1 to automatically refit when a new optimum is found
% % % opt_likreg.skipprof = 0; % skip profiling step; use boundaries from saved likreg set (1) or profiling (2)
% % % 
% % % par_better = calc_likregion(par_out,500,opt_likreg,opt_prof,opt_optim); 
% % % % Second entry is the number of accepted parameter sets to aim for. Use -1
% % % % here to use a saved set.
% % % 
% % % if isstruct(par_better) % if the profiling found a better optimum ...
% % %     print_par(par_better) % display on screen in formatted manner, so it can be copied into the code
% % %     calc_and_plot(par_better,opt_plot); % calculate model lines and plot them
% % %     return % stop here, and don't go into plotting with CIs
% % % end
% % % 
% % % %% Plot results with confidence intervals
% % % % The following code can be used to make a standard plot (the same as for
% % % % the fits), but with confidence intervals. Options for confidence bounds
% % % % on model curves can be set using opt_conf (see prelim_checks).
% % % % 
% % % % Use opt_conf.type to tell calc_conf which sample to use: 
% % % % -1) Skips CIs (zero does the same, and an empty opt_conf as well).
% % % % 1) Bayesian MCMC sample (default); CI on predictions are 95% ranges on 
% % % % the model curves from the sample 
% % % % 2) parameter sets from a joint likelihood region using the shooting 
% % % % method (limited sets can be used), which will yield (asymptotically) 95% 
% % % % CIs on predictions
% % % % 3) as option 2, but using the parameter-space explorer
% % % 
% % % opt_conf.type    = 2; % make intervals from 1) slice sampler, 2) likelihood region shooting, 3) parspace explorer
% % % opt_conf.lim_set = 2; % use limited set of n_lim points (1) or outer hull (2, likelihood methods only) to create CIs
% % % opt_conf.sens    = 0; % type of analysis 0) no sensitivities 1) corr. with state, 2) corr. with state/control, 3) corr. with relative change of state over time
% % % 
% % % out_conf = calc_conf(par_out,opt_conf);   % calculate confidence intervals on model curves
% % % calc_and_plot(par_out,opt_plot,out_conf); % call the plotting routine again to plot fits with CIs
% % % 
% % % % Here, we can also use the new plotting function for TKTD models. Even
% % % % though this is not a TKTD model, we can still plot the internal
% % % % concentration, with the treatments in separate panels.
% % % glo.locC       = [1 2]; % tell plot_tktd that our first and second state variable are internal concentrations to plot
% % % opt_tktd.repls = 0;     % plot individual replicates (1) or means (0)
% % % opt_tktd.min   = 1;     % set to 1 to show a dotted line for the control (lowest) treatment
% % % 
% % % plot_tktd(par_out,opt_tktd,opt_conf);
% % % 
% % % % %% Calculate LCx versus time
% % % % % Here, the LCx (by default the LC50) is calculated at several time points.
% % % % % LCx values are also printed on screen. If a sample from parameter space
% % % % % is available (e.g., from the slice sampler or the likelihood region), it
% % % % % can be used to calculate confidence bounds. 
% % % % % 
% % % % % Options for LCx (with confidence bounds) can be set using opt_ecx (see
% % % % % prelim_checks). Note that opt_conf.type=-1 skips CIs.
% % % % 
% % % % opt_conf.type    = 2; % make intervals from 1) slice sampler, 2)likelihood region, 3) parspace explorer
% % % % opt_conf.lim_set = 2; % for lik-region sample: use limited set of n_lim points (1) or outer hull (2) to create CIs
% % % % 
% % % % opt_ecx.id_sel    = [0   0 1]; % scenario to use from X0mat, scenario identifier, flag for ECx to use scenarios rather than concentrations
% % % % % opt_ecx.id_sel    = [6   6 1]; % scenario to use from X0mat, scenario identifier, flag for ECx to use scenarios rather than concentrations
% % % % % opt_ecx.id_sel    = [12 12 1]; % scenario to use from X0mat, scenario identifier, flag for ECx to use scenarios rather than concentrations
% % % % 
% % % % % This is the general method as the fast methods won't work for the full
% % % % % model (or at least: won't be faster than the general method).
% % % % opt_ecx.Feff      = [0.50]; % effect levels (>0 en <1), x/100 in ECx
% % % % opt_ecx.notitle   = 1; % set to 1 to suppress titles above ECx plots
% % % % Tend = [1:28]; % times at which to calculate LCx, relative to control
% % % % 
% % % % calc_ecx([],Tend,opt_ecx,opt_conf); % general method for ECx values