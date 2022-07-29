%% BYOM function derivatives.m (the model in ODEs)
%
%  Syntax: dX = derivatives(t,X,par,c,glo)
%
% This function calculates the derivatives for the full GUTS model system.
% Note that the survival probability dues to chemical stress is all
% calculated in <call_deri.html call_deri.m>. As input, it gets:
%
% * _t_   is the time point, provided by the ODE solver
% * _X_   is a vector with the previous value of the states
% * _par_ is the parameter structure
% * _c_   is the external concentration (or scenario number)
% * _glo_ is the structure with information (normally global)
%
% Time _t_ and scenario name _c_ are handed over as single numbers by
% <call_deri.html call_deri.m> (you do not have to use them in this
% function). Output _dX_ (as vector) provides the differentials for each
% state at _t_.
%
% * Author: Tjalling Jager
% * Date: September 2020
% * Web support: <http://www.debtox.info/byom.html>
% * Back to index <walkthrough_guts.html>

%  Copyright (c) 2012-2020, Tjalling Jager, all rights reserved.
%  This source code is licensed under the MIT-style license found in the
%  LICENSE.txt file in the root directory of BYOM. 

%% Start

function dX = derivatives(t,X,par,c,glo)

%% Unpack states
% The state variables enter this function in the vector _X_. Here, we give
% them a more handy name.

S  = X(glo.locS); % state is the survival probability at previous time point
Ci = X(glo.locC); % state is the internal concentration at previous time point
Di = X(glo.locD); % state is the scaled damage at previous time point

%% Unpack parameters
% The parameters enter this function in the structure _par_. The names in
% the structure are the same as those defined in the byom script file.
% The 1 between parentheses is needed as each parameter has 5 associated
% values.

ku   = par.ku(1);   % uptake rate constant, L * kg-1 d-1
ke   = par.ke(1);   % elimination rate constant
km   = par.km(1);     % formation rate of the metabolite
%kem  = par.kem(1);    % elimination rate of the metabolite

kr   = par.kr(1);   % damage repair rate constant
% mi   = par.mi(1);   % median of threshold distribution (used in call_deri)
% bi   = par.bi(1);   % killing rate (used in call_deri)
% Fs   = par.Fs(1);   % fraction spread of threshold distribution, (-) (used in call_deri)
hb   = par.hb(1);   % background hazard rate
T_A_tk  = par.T_A_tk(1);  % Arrhenius temperature for TK parameter, Kelvin
T_A_td  = par.T_A_td(1);  % Arrhenius temperature for TD parameter, Kelvin

% Unpack the experimental temperature based on the scenario
ref_temp = 293.15; % reference temperature in Kelvin (20 degrees celsius)
exp_temp = ref_temp; % initialise to pass error exp_temp (but it will be changed below)

% Define experimental temperatures for all scenarios 
exp_temp = glo.Temp_scen(2,c == glo.Temp_scen(1,:)); 

% Temperature correction for rates
ku_T = ku * exp( (T_A_tk / ref_temp) - (T_A_tk / exp_temp) );
ke_T = ke * exp( (T_A_tk / ref_temp) - (T_A_tk / exp_temp) );
km_T  = km * exp( (T_A_tk / ref_temp) - (T_A_tk / exp_temp) );
%kem_T = kem * exp( (T_A_tk / ref_temp) - (T_A_tk / exp_temp) );

kr_T = kr * exp( (T_A_td / ref_temp) - (T_A_td / exp_temp) );
hb_T = hb * exp( (T_A_td / ref_temp) - (T_A_td / exp_temp) );
%% Extract correct exposure for THIS time point
% Allow for external concentrations to change over time, either
% continuously, or in steps, or as a static renewal with first-order
% disappearance. For constant exposure, the code in this section is skipped
% (and could also be removed).

if glo.timevar(1) == 1 % if we are warned that we have a time-varying concentration ...
    c = read_scen(-1,c,t,glo); % use read_scen to derive actual exposure concentration
    % Note: this has glo as input to the function to save time!
    % the -1 lets read_scen know we are calling from derivatives (so need one c)   
end
%disp(c);
%% Calculate the derivatives
% This is the actual model, specified as a system of ODEs.

% dCi = ke * (Kiw * c - Ci); % first order bioconcentration
dCi = (ku_T * c - ke_T * Ci) - km_T * Ci ; % AMD: alternative writing of quation of first order bioconcentration
%dCm = km_T * Ci - kem_T * Cm;   % first-order metabolism

if glo.fastrep == 1 % is we assume that damage repair is infinitely fast
    dDi = dCi; % same change in damage as in internal conc.
else
    dDi = kr_T * (Ci - Di); % first order damage build-up from Ci (scaled)
end

dS = -hb_T * S; % only background hazard rate
% mortality due to the chemical is included in call_deri!

dX = zeros(size(X)); % initialise with zeros in right format
dX([glo.locS glo.locC glo.locD]) = [dS;dCi;dDi]; % collect derivatives in one vector


