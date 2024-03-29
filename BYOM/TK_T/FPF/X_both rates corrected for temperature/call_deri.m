%% BYOM function call_deri.m (calculates the model output)
%
%  Syntax: [Xout TE] = call_deri(t,par,X0v)
%
% This function calls the ODE solver to solve the system of differential
% equations specified in <derivatives.html derivatives.m>, or the explicit
% function(s) in <simplefun.html simplefun.m>. As input, it gets:
%
% * _t_   the time vector
% * _par_ the parameter structure
% * _X0v_   a vector with initial states and one concentration (scenario number)
%
% The output _Xout_ provides a matrix with time in rows, and states in
% columns. This function calls <derivatives.html derivatives.m>. The
% optional output _TE_ is the time at which an event takes place (specified
% using the events function). The events function is set up to catch
% discontinuities. It should be specified according to the problem you are
% simulating. If you want to use parameters that are (or influence) initial
% states, they have to be included in this function.
%
% * Author: Tjalling Jager
% * Date: January 2017
% * Web support: <http://www.debtox.info/byom.html>
% * Back to index <walkthrough_byom.html>

%% Start

function [Xout,TE,Xout2,zvd] = call_deri(t,par,X0v,glo)

% These outputs need to be defined, even if they are not used
Xout2    = []; % additional uni-variate output, predefine as empty
zvd      = []; % additional zero-variate output, predefine as empty

%% Initial settings
% This part extracts optional settings for the ODE solver that can be set
% in the main script (defaults are set in prelim_checks). The useode option
% decides whether to calculate the model results using the ODEs in
% <derivatives.html derivatives.m>, or the analytical solution in
% <simplefun.html simplefun.m>. Using eventson=1 turns on the events
% handling. Also modify the sub-function at the bottom of this function!
% Further in this section, initial values can be determined by a parameter
% (overwrite parts of X0), and zero-variate data can be calculated. See the
% example BYOM files for more information.

useode   = glo.useode; % calculate model using ODE solver (1) or analytical solution (2)
eventson = glo.eventson; % events function on (1) or off (0)
stiff    = glo.stiff; % set to 1 to use a stiff solver instead of the standard one

% Unpack the vector X0v, which is X0mat for one scenario
X0 = X0v(2:end); % these are the intitial states for a scenario

% if needed, calculate model values for zero-variate data from parameter
% set; these lines can be removed if no zero-variate data are used
if ~isempty(glo.zvd) % if there are zero-variate data defined (see byom_bioconc_extra)
    zvd       = glo.zvd; % copy zero-variate data structure to zvd
    zvd.BCF(3) = par.ku(1) / par.ke(1); % add model prediction as third value in zvd
else % if there are no zero-variate data defined (as in byom_bioconc_start)
    zvd       = []; % additional zero-variate output, output defined as empty matrix
end

%% Calculations
% This part calls the ODE solver (or the explicit model in <simplefun.html
% simplefun.m>) to calculate the output (the value of the state variables
% over time). There is generally no need to modify this part. The solver
% ode45 generally works well. For stiff problems, the solver might become
% very slow; you can try ode15s instead.

c  = X0v(1);     % the concentration (or scenario number)
t  = t(:);       % force t to be a row vector (needed when useode=0)

% specify options for the ODE solver
options = odeset; % start with default options
if eventson == 1
    options = odeset(options, 'Events',@eventsfun); % add an events function
end
options = odeset(options, 'RelTol',1e-4,'AbsTol',1e-7); % specify tightened tolerances
% options = odeset(options,'InitialStep',max(t)/1000,'MaxStep',max(t)/100); % specify smaller stepsize

TE = 0; % dummy for time of events

if useode == 1 % use the ODE solver to calculate the solution
    % call the ODE solver (use ode15s for stiff problems, escpecially for pulses)
    if isempty(options.Events) % if no events function is specified ...
        switch stiff
            case 0
                [~,Xout] = ode45(@derivatives,t,X0,options,par,c,glo);
            case 1
                [~,Xout] = ode113(@derivatives,t,X0,options,par,c,glo);
            case 2
                [~,Xout] = ode15s(@derivatives,t,X0,options,par,c,glo);
        end
    else % with an events functions ... additional output arguments for events:
        % TE catches the time of an event, YE the states at the event, and IE the number of the event
        switch stiff
            case 0
                [~,Xout,TE,YE,IE] = ode45(@derivatives,t,X0,options,par,c,glo);
            case 1
                [~,Xout,TE,YE,IE] = ode113(@derivatives,t,X0,options,par,c,glo);
            case 2
                [~,Xout,TE,YE,IE] = ode15s(@derivatives,t,X0,options,par,c,glo);
        end
    end
else % alternatively, use an explicit function provided in simplefun!
    Xout = simplefun(t,X0,par,c);
end

if isempty(TE) || all(TE == 0) % if there is no event caught
    TE = +inf; % return infinity
end

%% Events function
% Modify this part of the code if _eventson_=1 and _useode_=1. 
% This subfunction catches the 'events': in this case, it looks for the
% time where depuration starts.
%
% Note that the eventsfun has the same inputs, in the same sequence, as
% <derivatives.html derivatives.m>.

function [value,isterminal,direction] = eventsfun(t,X,par,c,glo)

nevents = 3;         % number of events that we try to catch

value       = zeros(nevents,1); % initialise with zeros
%if c == 17.27
%    value(1)   = t - 1.95; % catch depuration phase
%elseif c == 18
%    value(2)   = t - 1.95; % catch depuration phase
%elseif c == 19
%    value(3)   = t - 1.95; % catch depuration phase
%end
%value(1)    = X(1) - Ct;        % thing to follow is external concentration (state 1) minus threshold
isterminal  = zeros(nevents,1); % do NOT stop the solver at an event
direction   = zeros(nevents,1); % catch ALL zero crossing when function is increasing or decreasing
