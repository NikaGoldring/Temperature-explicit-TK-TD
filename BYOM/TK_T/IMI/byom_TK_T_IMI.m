%% BYOM, byom_calanus_2016_onecomp.m 
%
% * Author: Tjalling Jager 
% * Date: August 2018
% * Web support: <http://www.debtox.info/byom.html>
%
% BYOM is a General framework for simulating model systems in terms of
% ordinary differential equations (ODEs). The model itself needs to be
% specified in <derivatives.html derivatives.m>, and <call_deri.html
% call_deri.m> may need to be modified to the particular problem as well.
% The files in the engine directory are needed for fitting and plotting.
% Results are shown on screen but also saved to a log file (results.out).
%
% *The model:* An organism is exposed to a chemical in its surrounding
% medium. The animal accumulates the chemical according to standard
% one-compartment first-order kinetics.  
%
% *This script:* One-compartment TK of C2-naphthalene in _Calanus finmarchicus_.

%% MDA Changes: 
% Date: 10.04.2021
%
% Model extended to use datasets at different temperatures (see derivatives).
%
% Date: 21.07.2022
% Model changed to be used in BYOMv.6.2
% Model extended to account for biotransformation
%
%% Initial things
% Make sure that this script is in a directory somewhere *below* the BYOM
% folder.

clear, clear global % clear the workspace and globals
global DATA W X0mat % make the data set and initial states global variables 
global glo          % allow for global parameters in structure glo
diary off           % turn of the diary function (if it is accidentaly on)
set(0,'DefaultFigureWindowStyle','docked'); % collect all figure into one window with tab controls
%set(0,'DefaultFigureWindowStyle','normal'); % separate figure windows

pathdefine % set path to the BYOM/engine directory
glo.basenm  = mfilename; % remember the filename for THIS file for the plots
glo.saveplt = 0; % save all plots as (1) Matlab figures, (2) JPEG file or (3) PDF (see all_options.txt)

%% The data set
% Data are entered in matrix form, time in rows, scenarios (exposure
% concentrations) in columns. First column are the exposure times, first
% row are the concentrations or scenario numbers. The number in the top
% left of the matrix indicates how to calculate the likelihood:
%
% * -1 for multinomial likelihood (for survival data)
% * 0  for log-transform the data, then normal likelihood
% * 0.5 for square-root transform the data, then normal likelihood
% * 1  for no transformation of the data, then normal likelihood

% 7 degrees dataset
DATA{1,1} = [1	 17.71	17.71	17.71
           0.26  11.16	 6.68	19.29 % internal concentrations
           1.00	 32.56	44.23	28.37
           2.01	 64.18	64.22	57.20
           3.00	 54.03	43.32	59.00
           3.99	 46.85	70.16	50.28
           5.02	 57.44	50.88	62.34];
%DATA{1,2} = [];

% 18 degrees dataset   
% parent compound
DATA{2,1} = [1	 17.58	17.58	17.58
           0.27	 12.96	15.35	22.09
           0.99	 51.81	47.85	63.17
           1.95	 88.79	86.77	90.76
           2.88	 75.66	72.46	82.08
           3.87	 76.14	69.32	71.35
           4.99	 62.34	62.90	67.42];
% metabolite IMI-ole
DATA{2,2} = [1	    17.58	17.58	17.58
            0.27 	nan     nan     nan
            0.99	nan	    nan     nan
            1.95    nan     nan     nan 
            2.88	2.62	1.71	2.38
            3.87	1.87	2.14	2.13
            4.99	2.78	1.67	2.98];
        
% 24 degrees dataset
% parent compound
DATA{3,1} = [1	 17.57	17.57	17.57
           0.24	 16.95	23.40	23.83
           1.01	 69.93	74.25	64.62
           1.99	 93.05	95.68  112.35
           2.94	 68.56	92.55	68.71
           3.93	 69.47	65.31	67.15
           5.08	 55.05	65.18	57.53];      
% metabolite IMI-ole
DATA{3,2} = [ 1	    17.57	17.57	17.57
            0.24	nan     nan     nan		
            1.01	nan     nan     nan		
            1.99	2.96	3.29	3.59
            2.94	3.04	6.06	5.03
            3.93	4.40	4.28	4.01
            5.08	3.01	4.59	4.21 ];   

%% Initial values for the state variables
% Initial states, scenarios in columns, states in rows. First row are the
% 'names' of all scenarios.

X0mat = [17.71	17.58	17.57    % the scenarios (here nominal concentrations) 
          0      0       0       % initial values state 1 (internal concentrations)
          0      0       0];     % initial values state 2 (internal concentrations of metabolite)

glo.dep_time = [2.01	1.95	1.99]; % depuration times in days

glo.exp_temp = [280.15 291.15 297.15]; % experimental temperatures in Kelvin

glo.ref_temp = 293.15; % reference temperature in Kelvin (20 degrees celsius)

%% Initial values for the model parameters
% Model parameters are part of a 'structure' for easy reference. 

% syntax: par.name = [startvalue fit(0/1) minval maxval optional:log/normal scale (0/1)];
par.ke    = [0.1195    1 0.01 10 0];  % elimination rate constant, d-1
par.ku    = [3.177     1 0.01 10 0];  % uptake rate constant, L/kg/d
par.km    = [0.02191  1 1e-4 10 0];  % formation rate of the metabolite
par.kem   = [0.5188   1 0.01 10 0];  % elimination rate of the metabolite
par.T_A   = [2856  1 0.01 2e4 1];  % Arrhenius temperature, Kelvin

%% Zero-variate data and priors for Bayesian analyses
% Optionally, zero-variate data can be included. The corresponding model
% value needs to be calculated in call_deri.m, so modify that file too! You
% can make up your own parameter names in the structure _zvd_. Optionally,
% prior distributions can be specified for parameters, see the file
% calc_prior.m for the definition of the distributions. You must use the
% exact same names for the prior parameters as used in the _par_ structure.
% If you do not specify _zvd_ and/or _pri_ in your scripts, these options
% are simply not used in the analysis.
% 
% Note: the priors are also used for 'frequentist' analysis. They are
% treated as independent additional likelihood functions.
 
zvd.BCF = [20 0.5]; % example zero-variate data point for bioconcentration factor, L/kg, with normal s.d.
 
% % First element in pri is the choice of distribution.
% pri.Piw   = [2 116 121 118]; % triangular with min, max and center
% pri.ke    = [3 0.2 0.1];     % normal with mean and sd
% % If no prior is defined, the min-max bounds will define a uniform one.
% % Note that the prior is always defined on normal scale, also when the
% % parameter will be fitted on log scale.

%% Time vector and labels for plots
% Specify what to plot. If time vector glo.t is not specified, a default is
% used, based on the data set

% specify the y-axis labels for each state variable
glo.ylab{1} = ['internal concentration (',char(181),'g/kg)'];
glo.ylab{2} = ['internal metabolite concentration (',char(181),'g/kg)'];
% specify the x-axis label (same for all states)
glo.xlab    = 'time (day)';
glo.leglab1 = ''; % legend label before the 'scenario' number
glo.leglab2 = [char(181),'L']; % legend label after the 'scenario' number

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


% % %% Calculations and plotting
% % % Here, the functions are called that will do the calculation and the
% % % plotting. Note that calc_plot can provide all of the plotting information
% % % as output, so you can also make your own customised plots. This section,
% % % by default, makes a multiplot with all state variables (each in its own
% % % panel of the multiplot). When the model is fitted to data, output is
% % % provided on the screen (and added to the log-file results.out). The
% % % zero-variate data point is also plotted with its prediction (although the
% % % legend obscures it here).
% % %
% % % Options for the plotting can be set using opt_plot (see prelim_checks.m).
% % % Options for the optimisation routine can be set using opt_optim. Options
% % % for the ODE solver are part of the global glo.
% % %
% % % You can turn on the events function there too, to smoothly catch the
% % % discontinuity in the model. For the demo, the iterations were turned off
% % % (opt_optim.it = 0).
% % 
% % glo.eventson   = 1; % events function on (1) or off (0)
% % glo.useode     = 1; % calculate model using ODE solver (1) or analytical solution (0)
% % opt_optim.it   = 0; % show iterations of the optimisation (1, default) or not (0)
% % opt_plot.annot = 2; % extra subplot in multiplot for fits: 1) box with parameter estimates, 2) overall legend
% % opt_plot.bw    = 1; % if set to 1, plots in black and white with different plot symbols
% % 
% % % glo.diary = 'test_results.out'; % use a different name for the diary
% % 
% % % optimise and plot (fitted parameters in par_out)
% % par_out = calc_optim(par,opt_optim); % start the optimisation
% % calc_and_plot(par_out,opt_plot);     % calculate model lines and plot them
% % 
% % % save_plot(gcf,'fit_example',[],3) % save active figure as PDF (see save_plot for more options)
% % 
% % % return % stop here, and run analyses below manually later
% % 
% % %% Local sensitivity analysis
% % % Local sensitivity analysis of the model. All model parameters are
% % % increased one-by-one by a small percentage. The sensitivity score is
% % % by default scaled (dX/X p/dp) or alternatively absolute (dX p/dp).
% % %
% % % Options for the sensitivity can be set using opt_sense (see
% % % prelim_checks.m).
% %  
% % % % UNCOMMENT LINE(S) TO CALCULATE
% % % opt_sens.state = 2; % vector with state variables for the sensitivity analysis (0 does all, if sens>0)
% % % calc_localsens(par_out,opt_sens)
% % 
% % %% Profiling the likelihood
% % % By profiling you make robust confidence intervals for one or more of your
% % % parameters. Use the names of the parameters as they occurs in your
% % % parameter structure _par_ above. This can be a single string (e.g.,
% % % 'kd'), a cell array of strings (e.g., {'kd','ke'}), or 'all' to profile
% % % all fitted parameters. This example produces a profile for each parameter
% % % and provides the 95% confidence interval (on screen and indicated by the
% % % horizontal broken line in the plot).
% % %
% % % Notes: profiling for complex models is a slow process, so grab a coffee!
% % % If the profile finds a better solution, it breaks the analysis (as long
% % % as you keep the default opt_prof.brkprof=1) and displays the parameters
% % % for the new optimum on screen (and in results.out). For the NON-parallel
% % % version of calc_proflik, the new optimum is also immediately saves to the
% % % log file profiles_newopt.out. You can break of the analysis by pressing
% % % ctrl-c anytime, and use the values from the log file to restart (copy the
% % % better values into your _par_ structure). For parallel processing, saving
% % % would need more thought.
% % %
% % % Options for profiling can be set using opt_prof (see prelim_checks.m).
% % 
% % opt_prof.detail   = 1; % detailed (1) or a coarse (2) calculation
% % opt_prof.subopt   = 0; % number of sub-optimisations to perform to increase robustness
% % opt_prof.brkprof  = 2; % when a better optimum is located, stop (1) or automatically refit (2)
% % 
% % % % UNCOMMENT LINE(S) TO CALCULATE
% % % par_better = calc_proflik(par_out,'all',opt_prof,opt_optim);  % calculate a profile
% % % if ~isempty(par_better)                 % if the profiling found a better optimum ...
% % %     print_par(par_better) % display on screen in formatted manner, so it can be copied into the code
% % %     calc_and_plot(par_better,opt_plot); % calculate model lines and plot them
% % % end
% % 
% % %% Likelihood region
% % % Another way to make intervals on model predictions is to use a sample of
% % % parameter sets taken from the joint likelihood-based conf. region. This
% % % is done by the function calc_likregion.m. It first does profiling of all
% % % fitted parameters to find the edges of the region. Then, Latin-Hypercube
% % % shooting, keeping only those parameter combinations that are not rejected
% % % at the 95% level in a lik.-rat. test. The inner rim will be used for CIs
% % % on forward predictions.
% % %
% % % Options for the likelihood region can be set using opt_likreg (see
% % % prelim_checks.m). For the profiling part, use the options in opt_prof.
% % 
% % opt_prof.detail  = 1; % detailed (1) or a coarse (2) calculation
% % opt_prof.subopt  = 10; % number of sub-optimisations to perform to increase robustness
% % opt_prof.re_fit  = 1; % set to 1 to automatically refit when a new optimum is found
% % %opt_likreg.skipprof = 0; % skip profiling step; use boundaries from saved likreg set (1) or profiling (2)
% % 
% % par_better = calc_likregion(par_out,500,opt_likreg,opt_prof,opt_optim); 
% % % Second entry is the number of accepted parameter sets to aim for. Use -1
% % % here to use a saved set.
% % 
% % if isstruct(par_better) % if the profiling found a better optimum ...
% %     print_par(par_better) % display on screen in formatted manner, so it can be copied into the code
% %     calc_and_plot(par_better,opt_plot); % calculate model lines and plot them
% %     return % stop here, and don't go into plotting with CIs
% % end
% % 
% % %% Plot results with confidence intervals
% % % The following code can be used to make a standard plot (the same as for
% % % the fits), but with confidence intervals. Options for confidence bounds
% % % on model curves can be set using opt_conf (see prelim_checks).
% % % 
% % % Use opt_conf.type to tell calc_conf which sample to use: 
% % % -1) Skips CIs (zero does the same, and an empty opt_conf as well).
% % % 1) Bayesian MCMC sample (default); CI on predictions are 95% ranges on 
% % % the model curves from the sample 
% % % 2) parameter sets from a joint likelihood region using the shooting 
% % % method (limited sets can be used), which will yield (asymptotically) 95% 
% % % CIs on predictions
% % % 3) as option 2, but using the parameter-space explorer
% % 
% % opt_conf.type    = 2; % make intervals from 1) slice sampler, 2) likelihood region shooting, 3) parspace explorer
% % opt_conf.lim_set = 2; % use limited set of n_lim points (1) or outer hull (2, likelihood methods only) to create CIs
% % opt_conf.sens    = 0; % type of analysis 0) no sensitivities 1) corr. with state, 2) corr. with state/control, 3) corr. with relative change of state over time
% % 
% % out_conf = calc_conf(par_out,opt_conf);   % calculate confidence intervals on model curves
% % calc_and_plot(par_out,opt_plot,out_conf); % call the plotting routine again to plot fits with CIs
% % 
% % % % Here, we can also use the new plotting function for TKTD models. Even
% % % % though this is not a TKTD model, we can still plot the internal
% % % % concentration, with the treatments in separate panels.
% % % glo.locC       = [1 2]; % tell plot_tktd that our first and second state variable are internal concentrations to plot
% % % opt_tktd.repls = 0;     % plot individual replicates (1) or means (0)
% % % opt_tktd.min   = 0;     % set to 1 to show a dotted line for the control (lowest) treatment
% % % 
% % % plot_tktd(par_out,opt_tktd,opt_conf);