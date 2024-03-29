function [par_out,error_flag] = plot_grid_add(filename)

% Usage: [par_out,error_flag] = plot_grid_add(filename)
% 
% Function to combine MAT files for two chemicals (as generated by the
% parameter-space explorer) to make a MAT file for the mixture, assuming
% damage addition. This function makes a plot of mw x bw (SD) or Fs (IT)
% from multiple samples, to check visually whether they overlap (which
% would imply that damage addition is plausible). Two (or more) chemials
% can only be additive (in the strictest sense) if mw x bw is the same
% (SD), or Fs is the same (IT). Therefore, the parameter clouds for the
% individual chemicals are combined by only taking those combinations for
% which mw x bw or Fs is very close. This function only works properly if
% the MAT file was generated by the BYOM GUTS package for the standard
% reduced models (incl. ERA_special)! Output is also meant for the standard
% GUTS binary mixture package.
% 
% filename   cell array with strings for the name of MAT files to combine
% 
% Author     : Tjalling Jager
% Date       : December 2021
% Web support: <http://www.debtox.info/byom.html>

%  Copyright (c) 2012-2022, Tjalling Jager, all rights reserved.
%  This source code is licensed under the MIT-style license found in the
%  LICENSE.txt file in the root directory of BYOM. 

global glo2

SETTINGS_OPTIM = setup_settings(1); % load settings from file (use rough settings)
names = {'kd';'mw';'hb';'bw';'Fs'}; % names of the model parameters in the saved mat files
% this assumes that the standard GUTS package was used to create these files!

n_s = length(filename); % number of samples to compare
if n_s > 4 % for now, we only allow two files anyway, so this is for the future
    error_str = ('Cannot compare more than 4 samples at the moment');
    errordlg(error_str,'Error no. of sets'); % display a message box with the errors
    error_flag = 1; % signal main script that we ran into errors
    par_out    = [];
    return % simply stop as there is no useful input to work with
end

%% Run through sample to load files and see where 'good' sets are
% For damage addition, the parameters for chemical A and B are not
% independent. Under SD, the product of mw x bw must be the same, and under
% IT, Fs must be the same.

figure % create empty plot window for the comparison plot
hold on

% First, initialise matrices and cell arrays
sel   = nan(1,n_s);
BNDS  = cell(1,n_s);
COLL  = cell(1,n_s);
NAMES = cell(1,n_s);
HB    = cell(1,n_s);
MWBW  = cell(1,n_s);
FS    = cell(1,n_s);
plot_best = nan(1,n_s);
plot_neg  = nan(1,n_s);
plot_pos  = nan(1,n_s);

for i_s = 1:n_s % run through samples
    
    % Load sample <i_s>
    if exist([filename{i_s}],'file') == 2 % check if it exists first.
        load([filename{i_s}],'pmat','coll_all')
        % Load all of the saved information from the <calc_parspace> run.
    else % otherwise, produce an error (I don't think that is possible anymore)
        error(['There is no confidence set with filename ',filename{i_s},' in your working directory, so run calibration (with correct settings) first.'])
    end
    
    sel(i_s) = str2double(filename{i_s}(end-7)); % SD or IT (is part of the file name)
    
    % Extract useful information from <pmat>.
    ind_fit  = find(pmat(:,2) == 1); % indices to fitted parameters (vector)
    pmat_lim = pmat(ind_fit,:); % limit <pmat_lim> to fitted parameters
    ind_log  = find(pmat_lim(:,5) == 0); % indices to log-scale parameters in fitted <pmat> (vector)
    
    coll_all(:,ind_log) = 10.^coll_all(:,ind_log); % put back on normal scale, where needed
    names_lim           = names(ind_fit); % only keep names for fitted parameters
    
    % Chi2-criteria and indices to parameter sets to plot for inner and outer rim
    % We work here with the log-likelihood itself, so the chi2 criterion needs to be divided by 2
    chicrit_single = 0.5 * SETTINGS_OPTIM.crit_table(1,1);     % criterion for single-parameter CIs
    chicrit_prop   = 0.5 * SETTINGS_OPTIM.crit_prop(2);        % criterion for upper band of propagation band
    chicrit_prop2  = 0.5 * SETTINGS_OPTIM.crit_prop(1);        % criterion for lower band of propagation band
    
    % Find indices to sets within inner and outer rim
    ind_single = find(coll_all(:,end) < coll_all(1,end) + chicrit_single,1,'last'); % index to last element of sample that is still in inner rim
    
    % collect the results for further processing
    BNDS{i_s}  = pmat(ind_fit,[3 4]); % collect bounds used for original optimisation
    COLL{i_s}  = coll_all;
    NAMES{i_s} = names_lim;
    
    ind_hb  = find(strcmp(names,'hb')==1); % find location of <hb> in <pmat>
    HB{i_s} = pmat(ind_hb,1); % collect value for <hb>
    if pmat(ind_hb,2) == 1
        warning('off','backtrace')
        warning(['Background mortality was fitted in dataset ',filename{i_s},', this is ignored!'])
        warning('on','backtrace')
    end
    
    if sel(i_s) == 1 % for SD
        
        title('SD check for additivity')
        
        ind_mw = strcmp(names_lim,'mw'); % find location of <mw> in <pmat_lim>
        ind_bw = strcmp(names_lim,'bw'); % find location of <bw> in <pmat_lim>
        
        mwbw_all   = coll_all(:,ind_mw) .* coll_all(:,ind_bw); % product of mw and bw
        mwbw_inner = mwbw_all(1:ind_single,:); % points in the inner rim
        
        % plot comparison plot as best plus error bar
        plot_best(i_s) = mwbw_all(1); % best value
        plot_neg(i_s)  = plot_best(i_s) - min(mwbw_inner); % length of negative error bar
        plot_pos(i_s)  = max(mwbw_inner) - plot_best(i_s); % length of positive error bar        
        errorbar(i_s,plot_best(i_s),plot_neg(i_s),plot_pos(i_s),'ko','MarkerFaceColor','y','LineWidth',1)
        
%         % alternatively, plot points in the comparison plot
%         plot(i_s,mwbw_outer,'ko','MarkerFaceColor','w','MarkerEdgeColor','c') % plot all tries to continue with (or outer rim) as cyan open circles)
%         plot(i_s,mwbw_inner,'ko','MarkerFaceColor','g')
%         plot(i_s,mwbw_all(1),'ko','MarkerFaceColor','y')
        
        % collect the results for further processing
        MWBW{i_s} = mwbw_all;
        
    elseif sel(i_s) == 2 % for IT
        
        title('IT check for additivity')
        
        ind_Fs   = strcmp(names_lim,'Fs'); % find location of <Fs> in <pmat_lim>
        Fs_all   = coll_all(:,ind_Fs); % 
        Fs_inner = Fs_all(1:ind_single,:); % use the inner rim for those points (the ghost symbols)
        
        % plot comparison plot as best plus error bar
        plot_best(i_s) = Fs_all(1); % best value
        plot_neg(i_s)  = plot_best(i_s) - min(Fs_inner); % length of negative error bar
        plot_pos(i_s)  = max(Fs_inner) - plot_best(i_s); % length of positive error bar          
        errorbar(i_s,plot_best(i_s),plot_neg(i_s),plot_pos(i_s),'ko','MarkerFaceColor','y','LineWidth',1)
        
%         % alternatively, plot points in the comparison plot
%         plot(i_s,Fs_outer,'ko','MarkerFaceColor','w','MarkerEdgeColor','c') % plot all tries to continue with (or outer rim) as cyan open circles)
%         plot(i_s,Fs_inner,'ko','MarkerFaceColor','g')
%         plot(i_s,Fs_all(1),'ko','MarkerFaceColor','y')
        
        % collect the results for further processing
        FS{i_s} = Fs_all;
        
    end
    
    if sel(i_s)~=sel(1)
        error('all mat files for comparison need to be made with the same death mechanism')
    end
    
end

xlim([0 n_s+0.5])
xlabel('data set')
if sel == 1 % for SD
    ylabel('product mw x bw (1/d)')
elseif sel == 2 % for IT
    ylabel('Fs (-)')    
end
% modify tick on x-axis to show chemical A and B
xticks(1:n_s)
xticklabels({'A','B','C','D'}) % any extraneous entries will be ignored anyway
drawnow

% Test to see if the CIs overlap ... I think this works for binary
% mixtures, but not for tertiary and quaternary mixtures ...
plot_min = plot_best - plot_neg;
plot_max = plot_pos + plot_best;
ind_max = plot_max == max(plot_max); % finds highest upper CI of the comparison metric
% now the upper intervals of the other data sets must be above the lower
% interval of the max ...
ind_tst = plot_max ~= max(plot_max); % ones that do not have the highest upper CI
disp(' ')
if any(plot_max(ind_tst) < plot_min(ind_max))
    disp('It seems that the CIs do not overlap, indicating that damage addition is unlikely.')
else
    if n_s == 2
        disp('The CIs overlap, indicating that damage addition is a possibility.')
    else
        disp('Some CIs overlap, but this code cannot determine if they all do, so please check.')
    end
end

%% Create a new <coll_all> matrix for use with mixtures

limit      = 0.02; % factor within which the parameters must be to count as 'the same'
coll_all   = []; % start with empty <coll_all> matrix
glo2.names = {'hb';'kdA';'kdB';'mw';'bw';'Fs';'WB';'IAB'}; % names of the model parameters for new <coll_all> and <pmat>

% Start/check parallel pool
if glo2.n_cores > 0
    poolobj = gcp('nocreate'); % get info on current pool, but don't create one just yet
    if isempty(poolobj) % if there is no parallel pool ...
        parpool('local',glo2.n_cores) % create a local one with specified number of cores
    end
end

disp('Using parallel toolbox to combine MAT files. No progress will be shown.')

if sel(1) == 1 % for SD
    
    % find indices for all model parameters in names vector for each chemical
    ind_kdA = strcmp(NAMES{1},'kd');
    ind_kdB = strcmp(NAMES{2},'kd');
    ind_mwA = strcmp(NAMES{1},'mw');
    ind_bwA = strcmp(NAMES{1},'bw');
    ind_mwB = strcmp(NAMES{2},'mw');
    ind_bwB = strcmp(NAMES{2},'bw');
    
    coll_all_coll = cell(length(COLL{1}(:,1)),1);
    parfor iA = 1:length(COLL{1}(:,1)) % run through all elements for compound A
        
        % for each element in A, find where the product mw x bw is similar in B
        ind_ok = find(abs(MWBW{1}(iA)-MWBW{2})/MWBW{1}(iA) < limit);
        
        coll_add      = nan(length(ind_ok),6); % what we're adding to coll_all for set iA
        coll_add(:,1) = COLL{1}(iA,ind_kdA);     % collect kdA
        coll_add(:,2) = COLL{2}(ind_ok,ind_kdB); % collect kdB
        coll_add(:,3) = COLL{1}(iA,ind_mwA);     % collect mwA
        coll_add(:,4) = COLL{1}(iA,ind_bwA);     % collect bwA
        % weight factor will be the mean of the ratio mwA/mwB and bwB/bwA
        coll_add(:,5) = mean([COLL{1}(iA,ind_mwA)./COLL{2}(ind_ok,ind_mwB), COLL{2}(ind_ok,ind_bwB)./COLL{1}(iA,ind_bwA)],2);
        % min-log-likelihood is sum of that for the individual compounds
        coll_add(:,6) = COLL{1}(iA,end)+COLL{2}(ind_ok,end);
        
        coll_all_coll{iA} = coll_add;
        
    end
    
    for iA = 1:length(COLL{1}(:,1)) % run through all elements for compound A
        coll_all = cat(1,coll_all,coll_all_coll{iA}); % add to total coll_all
    end
    
elseif sel(1) == 2 % for IT

    % find indices for all model parameters in names vector for each chemical
    ind_kdA = strcmp(NAMES{1},'kd');
    ind_kdB = strcmp(NAMES{2},'kd');
    ind_mwA = strcmp(NAMES{1},'mw');
    ind_FsA = strcmp(NAMES{1},'Fs');
    ind_mwB = strcmp(NAMES{2},'mw');
    % ind_FsB = strcmp(NAMES{2},'Fs');
    
    coll_all_coll = cell(length(COLL{1}(:,1)),1);
    parfor iA = 1:length(COLL{1}(:,1)) % run through all elements for compound A
        ind_ok = find(abs(FS{1}(iA)-FS{2})/FS{1}(iA) < limit);
        coll_add = nan(length(ind_ok),6);
        coll_add(:,1) = COLL{1}(iA,ind_kdA);     % collect kdA
        coll_add(:,2) = COLL{2}(ind_ok,ind_kdB); % collect kdB
        coll_add(:,3) = COLL{1}(iA,ind_mwA);     % collect mwA
        coll_add(:,4) = COLL{1}(iA,ind_FsA);     % collect FsA
        % weight factor will be the ratio mwA/mwB
        coll_add(:,5) = COLL{1}(iA,ind_mwA)./COLL{2}(ind_ok,ind_mwB);
        % min-log-likelihood is sum of that for the individual compounds
        coll_add(:,6) = COLL{1}(iA,end)+COLL{2}(ind_ok,end);
        
        coll_all_coll{iA} = coll_add;
    end
    
    for iA = 1:length(COLL{1}(:,1)) % run through all elements for compound A
        coll_all = cat(1,coll_all,coll_all_coll{iA}); % add to total coll_all
    end
    
end

if isempty(coll_all)
    warning('No proper parameter combinations for addition can be found; no MAT file saved')
    error_flag = 1; % signal main script that we ran into errors
    par_out    = [];
    return % simply stop as there is no useful input to work with
end

coll_all      = sortrows(coll_all,size(coll_all,2));  % sort <coll_all> based on minloglik
chicrit_joint = 0.5 * SETTINGS_OPTIM.crit_table(5,1); % criterion for joint 95% CI or parameters
% Find indices to sets within inner and outer rim.
ind_fin95  = find(coll_all(:,end) < coll_all(1,end) + chicrit_joint,1,'last'); % index to last element of sample that is still in joint CI
ind_prop   = find(coll_all(:,end) < coll_all(1,end) + chicrit_prop,1,'last');  % index to last element of sample that is still in propagation band

coll_all(ind_fin95+1:end,:) = []; % remove the really bad ones

% If inner rim (incl. propagation set) is too large, we can downsample
if ind_prop > 6000
    ind_thin = randperm(ind_prop,ind_prop-5000); % downsample to keep an inner rim approx. 5000 sets (incl. propagation band)
    coll_all(1+ind_thin,:) = []; % remove some elements from outer rim! (but NOT the best value!)
    ind_fin95  = find(coll_all(:,end) < coll_all(1,end) + chicrit_joint,1,'last'); % index to last element of sample that is still in joint CI
    ind_prop   = find(coll_all(:,end) < coll_all(1,end) + chicrit_prop,1,'last');  % index to last element of sample that is still in propagation band
end

% we don't need so many in the outer rim, so we can downsample to decrease
% size of the mat file produced
thin_sz  = ind_fin95-ind_prop;
ind_thin = randperm(thin_sz,max(1,thin_sz-ind_prop*3)); % downsample to keep an outer rim approx. 3x the propagation set
coll_all(ind_prop+ind_thin,:) = []; % remove some elements from outer rim!

ind_fin95  = find(coll_all(:,end) < coll_all(1,end) + chicrit_joint,1,'last');  % index to last element of sample that is still in joint CI
ind_prop   = find(coll_all(:,end) < coll_all(1,end) + chicrit_prop,1,'last');   % index to last element of sample that is still in upper edge propagation band
ind_prop2  = find(coll_all(:,end) < coll_all(1,end) + chicrit_prop2,1,'last');  % index to last element of sample that is still in lower edge propagation band
ind_single = find(coll_all(:,end) < coll_all(1,end) + chicrit_single,1,'last'); % index to last element of sample that is still in inner rim

%% Create new <pmat> for mixture

pmat      = zeros(8,5); % create a new <pmat> for the mixture addition case
pmat(1,:) = [mean([HB{1} HB{2}])  0    0    1 1]; % for now, take <hb> as mean for two data sets!
if pmat(1,1) > 1e-6 % best to warn that this is happening
    warning('off','backtrace')
    warning('Background hazard hb is set to MEAN for two data sets!')
    warning('on','backtrace')
end

% fill <pmat> with information from new <coll_all> for parameters common to SD and IT
pmat(2,:) = [coll_all(1,1) 1 min(coll_all(:,1)) max(coll_all(:,1)) 0]; % kdA
pmat(3,:) = [coll_all(1,2) 1 min(coll_all(:,2)) max(coll_all(:,2)) 0]; % kdB
pmat(4,:) = [coll_all(1,3) 1 min(coll_all(:,3)) max(coll_all(:,3)) 0]; % mw
pmat(7,:) = [coll_all(1,5) 1 min(coll_all(:,5)) max(coll_all(:,5)) 0]; % WB
pmat(8,:) = [0     0 -100  100 1]; % IAB

% Make the bounds a bit wider, but not wider than the original bounds in the saved <pmat>
bnds_extra = 1.5; % what to divide/multiply min-max bounds with to make them a bit wider than sample
pmat([2 3],3) = max(pmat([2 3],3)/bnds_extra,[BNDS{1}(ind_kdA,1);BNDS{2}(ind_kdB,1)]); % min bounds for kdA and kdB
pmat([2 3],4) = min(pmat([2 3],4)*bnds_extra,[BNDS{1}(ind_kdA,2);BNDS{2}(ind_kdB,2)]); % max bounds for kdA and kdB
pmat(4,3) = max(pmat(4,3)/bnds_extra,BNDS{1}(ind_mwA,1)); % min bounds for mw
pmat(4,4) = min(pmat(4,4)*bnds_extra,BNDS{1}(ind_mwA,2)); % max bounds for mw
pmat(7,3) = pmat(7,3)/bnds_extra; % min bounds for WB
pmat(7,4) = pmat(7,4)*bnds_extra; % max bounds for WB

% fill <pmat> with information from new <coll_all> for parameters NOT common to SD and IT
if sel(1) == 1 % for SD
    pmat(5,:) = [coll_all(1,4) 1 min(coll_all(:,4)) max(coll_all(:,4)) 0]; % bw
    pmat(6,:) = [3     0    1  20 1]; % Fs
    pmat(5,3) = max(pmat(5,3)/bnds_extra,BNDS{1}(ind_bwA,1)); % min bounds for bw
    pmat(5,4) = min(pmat(5,4)*bnds_extra,BNDS{1}(ind_bwA,2)); % min bounds for bw
elseif sel(1) == 2
    pmat(5,:) = [1e3 0 1e-6 1e6 1]; % bw
    pmat(6,:) = [coll_all(1,4) 1 min(coll_all(:,4)) max(coll_all(:,4)) 0]; % Fs
    pmat(6,3) = max(pmat(6,3)/bnds_extra,BNDS{1}(ind_FsA,1)); % min bounds for Fs
    pmat(6,4) = min(pmat(6,4)*bnds_extra,BNDS{1}(ind_FsA,2)); % min bounds for Fs
end

ind_fit = find(pmat(:,2) == 1); % indices to fitted parameters (vector)

% bounds that are smaller than a factor 10 can be fitted on normal scale
chk_bnds = pmat(:,2) ==1 & (pmat(:,4)./pmat(:,3))<10; % which elements are smaller than 10?
pmat(chk_bnds,5) = 1; % put them to normal scale

%% Display, plot and final things

WRAP.glo  = [];
WRAP.glo2 = glo2;

% display par in a formatted way so they can be directly copied into the code of the main script
par_out = packunpack(2,[],pmat,WRAP);
disp(' ')
if sel(1) == 1
    disp('  Damage addition for stochastic death (SD)')
else
    disp('  Damage addition for individual tolerance (IT)')
end
disp('You can copy-paste following lines into script for mixture predictions')
disp('if you like to run the parameter-space explorer on the mixture data.')
print_par(par_out) % this prints out the optimised parameter values in a

pmat_lim = pmat(ind_fit,:); % limited parameter matrix for fitted parameters
ind_log  = find(pmat_lim(:,5) == 0); % indices to log-scale parameters in fitted <pmat> (vector)

% improvise a <pmat_print>, though without any CIs, so it can be saved
pmat_print      = zeros(length(ind_fit),8);
pmat_print(:,1) = pmat_lim(:,1); % copy best-fit values
% add CIs for the model parameters, estimated as the edges of the
% propagation set; this is a good approximation
pmat_print(:,2) = (min(coll_all(1:ind_prop,1:end-1),[],1))'; 
pmat_print(:,3) = (max(coll_all(1:ind_prop,1:end-1),[],1))';

disp(' ')
disp('=================================================================================')
disp('Results of the parameter-space combination for mixture predictions')
disp('=================================================================================')
if sel(1) == 1
    disp('Damage addition for stochastic death (SD)')
else
    disp('Damage addition for individual tolerance (IT)')
end
disp(['   Chemical A is taken from: ',filename{1}])
disp(['   Chemical B is taken from: ',filename{2}])
disp(['   Sample: ',num2str(ind_fin95),' sets in joint CI and ',num2str(ind_single),' in inner CI.'])
disp(['   Propagation set: ',num2str(1+ind_prop - ind_prop2),' sets will be used for error propagation.'])

FVAL = coll_all(1,end);
fprintf('   Minus log-likelihood has reached the value %#1.2f (AIC=%#1.2f). \n',FVAL,2*length(ind_fit)+2*FVAL)

disp('Approx. best estimates and 95% CIs on fitted parameters')
disp_pmat(pmat,pmat_print,WRAP); % call dedicated function for printing the <pmat>
disp('Confidence intervals are estimated from the sample (inner rim, incl. propagation band)')
disp('so they will generally slightly exaggerate the true intervals.')

% put <pmat> on log-scale for log-parameters (this is how BYOM uses it)
pmat(pmat(:,2)==1 & pmat(:,5)==0,1) = log10(pmat(pmat(:,2)==1 & pmat(:,5)==0,1));

% put <coll_all> on log-scale for log-parameters
coll_all(:,ind_log) = log10(coll_all(:,ind_log));

coll_prof_pruned = cell(1,length(ind_fit));
for i_p = 1:length(ind_fit) % trick to make sure profiles are plotted (but without refinement red line)
    coll_prof_pruned{i_p} = [NaN NaN];
end

% And make a plot of the final results as parameter-space plot
figh = plot_grid(pmat,[],coll_all,coll_prof_pruned,[],SETTINGS_OPTIM,WRAP); 

% user can give a filename for saving, but I try to provide a good default
filenm = ['ADD_',filename{1}(11:end-12),'+',filename{2}(11:end)];

% % Instead of the line above, feel free to uncomment the code below to use
% % the Matlab GUI for entering a filename for saving.
% filenm = uiputfile(['ADD_',filename{1}(11:end-12),'+',filename{2}(11:end)],'Save new project file for mixture of chemical A and B'); % use Matlab GUI to select name for saving MAT file
% if numel(filenm) == 1 && filenm(1) == 0 % if cancel is pressed ...
%     error_flag = 1; % signal main script that we ran into errors
%     return % simply stop as the user does not want to save
% end

saveplt = 0; % for now, don't save (since this function is not normally called from a standard byom script, the glo.saveplot will not be defined)
figure(figh)  % use existing handle and make sure that the multiplot is the current plot
% Add a title to the plot.
axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping','off');
text(0.5, 1,['File: ',filenm,' (',date,')'],'Interpreter','none','HorizontalAlignment','center','VerticalAlignment', 'top');
if saveplt > 0
    save_plot(figh,['parspace_addmixpred_',glo.basenm]) % save parspace plot in output folder
end
snapnow % for publishing Matlab scripts: make plot appear in the context of the screen output

% save new mat file
par   = par_out;
names = glo2.names ;
save(filenm,'par','pmat','coll_all','pmat_print','coll_prof_pruned','names')
