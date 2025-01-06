clear;
close all;
clc;

%% user inputs

% name and address of batch file
fadd_batch = 'inputs\tandem_params_new.xlsx';

% variables to be imported
varnames = {'dist_odias', 'dat'};

% on-demand function to correct bias in AAC classification
coef_rho = @(x) -3.598 * x.^(-0.8376) + 1.27;


%% read batch file

% set data type and import excel file
opts_read = detectImportOptions(fadd_batch);
opts_read = setvartype(opts_read, 'char');
intab = readtable(fadd_batch, 'UseExcel', true);

% sort out batch file parameters
n_dat = intab.n_dat(1); % total number of data given by the user
ii0 = 1 : n_dat; % indices of the input data
fdir = intab.fdir(ii0); % folder addresses of MATLAB workspaces that...
    % ...contain inverted data
fname = intab.fname(ii0); % filenames for MATLAB workspaces
clr = intab.clr(ii0); % assign colors to plot markers and lines
mrk = intab.mrk(ii0); % assign marker symbols
linstl = intab.linstl(ii0); % assign line styles for plots
test_date = intab.date(ii0); % date of measurement
test_condition = intab.condt(ii0); % condition of measurement
incld = intab.incld(ii0); % determine whether a dataset would be...
    % ...post-processed
ii = ii0(incld); % indices of data selected by the user for post-processing

%% initialize groups assigned in batch file (if requested)
if ismember('group', fieldnames(intab))

    ind0_grp = intab.group(1:n_dat); % read data grouping indices

    % find unique group indices
    ind_grp_uniq = unique(ind0_grp(ii));

    n_grp = numel(ind_grp_uniq); % number of unique groups
    
    % initialize placeholders for indices of datasets corresponding to...
        % ...each group
    ind_grp = cell(n_grp, 1);
    
    % find members of each group
    for k = 1 : n_grp
        ind_grp{k} = find((ind0_grp == ind_grp_uniq(k)));
        ind_grp{k} = ind_grp{k}(ismember(ind_grp{k},ii));
    end

    % initialize a strcuture for the group
    dist_grp = struct('d_mode', cell(n_grp, 1),...
    'd_gm', cell(n_grp, 1), 'rho_eff', cell(n_grp, 1));

end

% color for group plots
clr2 = {{'#8D493A', '#DC8686'}, {'#537188', '#7EACB5'}};
mrk2 = {'^', 'o'};

%% initialize universal correlation for rho_eff vs. dm

D_m = 2.48; % exponent
rho_eff_100 = 510; % pefactor
dm_lim_uc = [1e0 2e4];  % limits on the mobility diameter
n_dm_uc = 1e4; % number of data
uc = @(y) rho_eff_100 * (y / 100) .^ (D_m - 3); % on-demand function...
% ...for the forward correlation in the mass-mobility domain...
% ...(rho_eff in [kg/m3] as a function of dm in [nm])
r_uc = (dm_lim_uc(2) / dm_lim_uc(1)) ^ (1 / (n_dm_uc - 1));
dm_uc = dm_lim_uc(1) * ones(n_dm_uc,1) .* r_uc .^ (((1 : n_dm_uc) - 1)');
rho_eff_uc = uc(dm_uc);

%% load inverted 1d tandem distributions and calculate effective density

% allocate space to placeholders for ensemble of inverted size...
    % ...distributions to be loaded
dists = cell(n_dat,1);

% initialize the figure for effective densities segregated by day
f1 = figure(1);
f1.Position = [50, 50, 500, 600];
set(f1, 'color', 'white');

% initialize plot and legend placeholders
n_ii = length(ii);
plt1 = cell(n_dat + 1, 1);
lgdtxt1 = cell(n_dat + 1, 1);

% plot universal correlation
plt1{end} = plot(dm_uc, rho_eff_uc, 'Color', [0.4940 0.1840 0.5560],...
    'LineStyle', '-.', 'LineWidth', 2);
lgdtxt1{end} = 'Olfert \& Rogak (2019)';
hold on

for i = ii
    
    % import saved workspace data for each condition
    fadd_wsp = strcat(fdir{i}, '\', fname{i}, '.mat');
    load(fadd_wsp, varnames{:})
    
    % find the highest peak among the modes at each AAC setpoint
    dist_odias = hn.selectpeak(dist_odias);
    
    for j = 1 : size(dist_odias,1) 

        % calculate effective density for each setpoint
        dist_odias(j).rho_eff = DAT.RHO_EFF(dist_odias(j).d_mode,...
            dat(j).da) ./ coef_rho(dist_odias(j).d_mode);
        
        % store aerodynamic setpoint
        dist_odias(j).da = dat(j).da;
    
    end

    % store the processed data
    dists{i} = dist_odias;
    
    % plot scatters of effective density vs. mobility diameter at each day
    plt1{i} = scatter(cat(1,dists{i}.d_mode), cat(1,dists{i}.rho_eff),...
        15, hex2rgb(clr{i}), mrk{i}, 'LineWidth', 1.5);
    hold on
    
    % make a description for the dataset
    lgdtxt1{i} = strcat(test_date(i), ',', {' '}, test_condition(i));
    
    % clear redundant variables
    clear dist_odias dat rho_eff

end

% set apprearances
box on
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12,...
    'TickLength', [0.02 0.02], 'XScale', 'log', 'YScale', 'log')
xlabel('$d_\mathrm{m}$ [nm]', 'interpreter', 'latex',...
    'FontSize', 16)
ylabel('$\rho_\mathrm{eff} [\mathrm{kg}/\mathrm{m}^3]$',...
    'interpreter', 'latex', 'FontSize', 16)

%% present data in groups if requested by user

if exist('ind0_grp', 'var') && ~isempty(ind0_grp)

    iii = []; % initialize index for legends of figure 1

    % initialize the figure for ensemble groups
    f2 = figure(2);
    f2.Position = [150, 150, 500, 500];
    set(f2, 'color', 'white');
    
    % initialize plot and legend placeholders
    plt2 = cell(n_grp + 1, 1);
    lgdtxt2 = cell(n_grp + 1, 1);

    % plot universal correlation
    plt2{end} = plot(dm_uc, rho_eff_uc, 'Color', [0.4940 0.1840 0.5560],...
        'LineStyle', '-.', 'LineWidth', 2);
    lgdtxt2{end} = 'Olfert \& Rogak (2019)';
    hold on
    
    % loop over the groups
    for k = 1 : n_grp
        
        for kk = 1 : length(ind_grp{k})
            
            % compile modes, means and densities into different groups
            dist_grp(k).d_mode = [dist_grp(k).d_mode,...
                dists{ind_grp{k}(kk)}.d_mode];
            dist_grp(k).d_gm = [dist_grp(k).d_gm,...
                dists{ind_grp{k}(kk)}.d_gm];
            dist_grp(k).rho_eff = [dist_grp(k).rho_eff,...
                dists{ind_grp{k}(kk)}.rho_eff];
            
            % sort out indices for order of figure 1 legend
            iii = [iii, ind_grp{k}(kk)];
            
        end
              
        % plot the grouped effective density datasets
        plt2{k} = scatter(dist_grp(k).d_mode, dist_grp(k).rho_eff, 15,...
            hex2rgb(clr2{k}{1}), mrk2{k}, 'LineWidth', 1.5);
        hold on

        lgdtxt2{k} = strcat('Group', {' '}, num2str(k));

    end
    
    % set apprearances
    box on
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12,...
        'TickLength', [0.02 0.02], 'XScale', 'log', 'YScale', 'log')
    xlabel('$d_\mathrm{m}$ [nm]', 'interpreter', 'latex',...
        'FontSize', 16)
    ylabel('$\rho_\mathrm{eff} [\mathrm{kg}/\mathrm{m}^3]$',...
        'interpreter', 'latex', 'FontSize', 16)
    xlim([0.8 * min(cat(2,dist_grp.d_mode))...
        1.2 * max(cat(2,dist_grp.d_mode))])
    ylim([0.8 * min(cat(2,dist_grp.rho_eff))...
        1.2 * max(cat(2,dist_grp.rho_eff))])
    legend(cat(2, plt2{:}), cat(2, lgdtxt2{:}), 'interpreter', 'latex',...
    'FontSize', 12, 'Location', 'southwest');
    
    % adjust the bounds in non-grouped figure
    figure(f1)
    xlim([0.8 * min(cat(2,dist_grp.d_mode))...
        1.2 * max(cat(2,dist_grp.d_mode))])
    ylim([0.8 * min(cat(2,dist_grp.rho_eff))...
        1.2 * max(cat(2,dist_grp.rho_eff))])

    % sort out order of legends in non-grouped figure
    legend(cat(2, plt1{iii}, plt1{end}), cat(2, lgdtxt1{iii}, lgdtxt1{end}),...
    'interpreter', 'latex', 'FontSize', 10, 'Location', 'southoutside',...
    'NumColumns', 2);

else

    % just print legend for the non-grouped figure
    legend(cat(2, plt1{:}), cat(2, lgdtxt1{:}), 'interpreter', 'latex',...
    'FontSize', 10, 'Location', 'southoutside', 'NumColumns', 2);

end



