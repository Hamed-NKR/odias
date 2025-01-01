clear;
close all;
clc;

%% user inputs

% name and address of batch file
fadd_batch = 'inputs\tandem_params_new.xlsx';

% variables to be imported
varnames = {'dist_odias', 'dat'};

%% read batch file

% set data type and import excel file
opts_read = detectImportOptions(fadd_batch);
opts_read = setvartype(opts_read, 'char');
intab = readtable(fadd_batch, 'UseExcel', true);

% sort out batch file parameters
n_dat = intab.n_dat(1); % total number of data given by the user
ii = 1 : n_dat; % indices of the input data
incld = intab.incld(ii); % determine whether a dataset would be...
    % ...post-processed
ii = ii(incld); % indices of data selected by the user for post-processing
fdir = intab.fdir(ii); % folder addresses of MATLAB workspaces that...
    % ...contain inverted data
fname = intab.fname(ii); % filenames for MATLAB workspaces
clr = intab.clr(ii); % assign colors to plot markers and lines
mrk = intab.mrk(ii); % assign marker symbols
linstl = intab.linstl(ii); % assign line styles for plots

% initialize groups assigned in batch file
if ismember('group', fieldnames(intab))

    ind0_grp = intab.group(ii); % read data grouping indices

    % find unique group indices
    ind_grp_uniq = unique(ind0_grp(ii));

    n_grp = numel(ind_grp_uniq); % number of unique groups
    
    % initialize placeholders for indices of datasets corresponding to...
        % ...each group
    ind_grp = cell(n_grp, 1);
    
    % find members of each group
    for k = 1 : n_grp
        ind_grp{k} = find(ind0_grp(ii) == ind_grp_uniq(k));
    end

    % initialize a strcuture for the group
    dist_grp = struct('d_mode', cell(n_grp, 1),...
    'd_gm', cell(n_grp, 1), 'rho_eff', cell(n_grp, 1));

end

% color for group plots
clr2 = {{'#8D493A', '#DC8686'}, {'#537188', '#7EACB5'}};

%% load inverted 1d tandem distributions and calculate effective density

% allocate space to placeholders for ensemble of inverted size...
    % ...distributions to be loaded
dists = cell(n_dat,1);

for i = ii
    
    % import saved workspace data for each condition
    fadd_wsp = strcat(fdir{i}, '\', fname{i}, '.mat');
    load(fadd_wsp, varnames{:})
    
    % find the highest peak among the modes at each AAC setpoint
    dist_odias = hn.selectpeak(dist_odias);
    
    % calculate effective density for each setpoint
    for j = 1 : size(dist_odias,1) 
        dist_odias(j).rho_eff = DAT.RHO_EFF(dist_odias(j).d_mode,...
            dat(j).da);
        dist_odias(j).da = dat(j).da;
    end

    % store the processed data
    dists{i} = dist_odias;
    
    % clear redundant variables
    clear dist_odias dat rho_eff

end

%% present data in groups if requested by user

if exist('ind0_grp', 'var') && ~isempty(ind0_grp)

    % initialize size distribution figure for groups
    f2 = figure(2);
    f2.Position = [100, 100, 500, 500];
    set(f2, 'color', 'white');
    
    % initialize plot and legend placeholders
    plt2 = cell(n_grp,1);
    lgdtxt2 = cell(n_grp,1);

    % loop over the groups
    for k = 1 : n_grp
        
        for kk = 1 : length(ind_grp{k})
            
            % compile modes, means and densities into different groups
            dist_grp(k).d_mode = [dist_grp(k).d_mode;...
                dists{ind_grp{k}(kk)}.d_mode];
            dist_grp(k).d_gm = [dist_grp(k).d_gm;...
                dists{ind_grp{k}(kk)}.d_gm];
            dist_grp(k).rho_eff = [dist_grp(k).rho_eff;...
                dists{ind_grp{k}(kk)}.rho_eff];
            
        end
              
        % plot the average distributions and error bounds
        plt2{k} = scatter(dist_grp(k).d_mode, dist_grp(k).rho_eff, 12,...
            hex2rgb(clr2{k}{1}), mrk{k}, 'LineWidth', 1);
        hold on

        lgdtxt2{k} = strcat('Group', {' '}, num2str(k));

    end
    
    % set apprearances
    box on
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12,...
        'TickLength', [0.02 0.02], 'XScale', 'log')
    xlabel('$d_\mathrm{m}$ [nm]', 'interpreter', 'latex',...
        'FontSize', 16)
    ylabel('$\rho_\mathrm{eff}) [\mathrm{kg}/\mathrm{m}^3]$',...
        'interpreter', 'latex', 'FontSize', 16)
    xlim([15,1000])
    legend(cat(2, plt2{:}), cat(2, lgdtxt2{:}), 'interpreter', 'latex',...
    'FontSize', 12, 'Location', 'northeast');
    
end


