%% BYOM function derivatives.m (the model in ODEs)
%
%  Syntax: dX = derivatives(t,X,par,c,glo)
%
% This function calculates the derivatives for the model system. It is
% linked to the script files <byom_bioconc_extra.html byom_bioconc_extra.m>
% and <byom_bioconc_start.html byom_bioconc_start.m>. As input,
% it gets:
%
% * _t_   is the time point, provided by the ODE solver
% * _X_   is a vector with the previous value of the states
% * _par_ is the parameter structure
% * _c_   is the external concentration (or scenario number)
% * _glo_ is the structure with information (normally global)
% 
% Note: _glo_ is now an input to this function rather than a global. This
% makes the code considerably faster.
%
% Time _t_ and scenario name _c_ are handed over as single numbers by
% <call_deri.html call_deri.m> (you do not have to use them in this
% function). Output _dX_ (as vector) provides the differentials for each
% state at _t_.
%
% * Author: Tjalling Jager
% * Date: May 2020
% * Web support: <http://www.debtox.info/byom.html>
% * Back to index <walkthrough_byom.html>
% 
%  Copyright (c) 2012-2020, Tjalling Jager, all rights reserved.
%  This source code is licensed under the MIT-style license found in the
%  LICENSE.txt file in the root directory of BYOM. 

%% Start
% Note that, in this example, variable _c_ (the scenario identifier) is not
% used in this function. The treatments differ in their initial exposure
% concentration that is set in X0mat. Also, the structure _glo_ is not used
% here.

function dX = derivatives(t,X,par,c,glo)

%% Unpack states
% The state variables enter this function in the vector _X_. Here, we give
% them a more handy name.

ref_temp = glo.ref_temp;
if c == 17.71 % ADJUST!
    exp_temp = glo.exp_temp(1);
    dep_time = glo.dep_time(1); 
end

if c == 17.58 % ADJUST!
    exp_temp = glo.exp_temp(2);
    dep_time = glo.dep_time(2); 
end

if c == 17.57 % ADJUST!
    exp_temp = glo.exp_temp(3);
    dep_time = glo.dep_time(3); 
end

if t< dep_time  %AMD: modify for respective depuration time!
    Cw = c; % AMD: c is the scenario identifyer, here this refers to the external concentration 
    
else 
    Cw = 0;
end

Ci = X(1); % internal concentration
Cm = X(2);     % state 2: metabolite concentration

%% Unpack parameters
% The parameters enter this function in the structure _par_. The names in
% the structure are the same as those defined in the byom script file.
% The 1 between parentheses is needed as each parameter has 5 associated
% values.

% AMD:
ku   = par.ku(1);     % uptake rate constant, L * kg-1 d-1
ke   = par.ke(1);     % elimination rate constant, d-1
% Piw  = par.ku(1)/ke;  % bioconcentration factor, L/kg
km   = par.km(1);     % formation rate of the metabolite
kem  = par.kem(1);    % elimination rate of the metabolite
T_A  = par.T_A(1);    % Arrhenius temperature, Kelvin

% Correct rates for temperature 
ku_T  = ku * exp( (T_A / ref_temp) - (T_A / exp_temp(1)) );
ke_T  = ke * exp( (T_A / ref_temp) - (T_A / exp_temp(1)) );
km_T  = km * exp( (T_A / ref_temp) - (T_A / exp_temp(1)) );
kem_T = kem * exp( (T_A / ref_temp) - (T_A / exp_temp(1)) );

% % Optionally, include the threshold concentration as a global (see script)
% Ct   = glo.Ct;      % threshold external concentration that stops degradation, mg/L

%% Calculate the derivatives
% This is the actual model, specified as a system of two ODEs:
%
% $$ \frac{d}{dt}C_w=-k_d C_w \quad \textrm{as long as } C_w>C_t $$
%
% $$ \frac{d}{dt}C_i=k_e(P_{iw}C_w-C_i) $$

%dCi = ke * (Piw * Cw - Ci); % first order bioconcentration
dCi = (ku_T * Cw - ke_T * Ci) - km_T * Ci ; % AMD: alternative writing of quation of first order bioconcentration
dCm = km_T * Ci - kem_T * Cm;   % first-order metabolism

dX = [dCi; dCm]; % collect derivatives in vector dX
