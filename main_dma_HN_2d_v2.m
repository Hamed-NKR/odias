clear;
close all;
clc;
warning('off')

%% user inputs

% name and address of batch file
fadd_batch = 'inputs\tandem_params_new.xlsx';

% variables to be imported
varnames = {'dist_odias', 'dat'};

% on-demand function to correct bias in AAC classification
coef_rho = @(x) -3.598 * x.^(-0.8376) + 1.27;

% polynomial degree fit
fittype = {'poly3', 'poly3', 'poly4', 'poly4'};

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
mrksz = intab.mrksz(ii0); % assign marker size
linstl = intab.linstl(ii0); % assign line styles for plots
test_date = intab.date(ii0); % date of measurement
test_condition0 = intab.condt(ii0); % condition of measurement
grp = intab.group(ii0);
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
    dist_grp = struct('d_mode', cell(n_grp, 1), 'da', cell(n_grp, 1),...
    'd_gm', cell(n_grp, 1), 'rho_eff', cell(n_grp, 1),...
    'sigma_g', cell(n_grp, 1));

end

% color for group plots
clr2 = {{'#8D493A', '#DC8686'}, {'#537188', '#7EACB5'},...
    {'#8174A0', '#A888B5'}, {'#659287', '#B1C29E'}};
mrk2 = {'^', 'o', 's', 'h'};
mrksz2 = [20, 20, 30, 25];

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
f1.Position = [50, 50, 800, 600];
set(f1, 'color', 'white');
t1 = tiledlayout(1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

% initialize placeholders for plots and legends 
n_ii = length(ii);
plt1 = cell(n_dat + 1, 1);
plt3 = cell(n_dat + 1, 1);
lgdtxt1 = cell(n_dat + 1, 1);

% plot the universal correlation
nexttile(1)
plt1{end} = plot(dm_uc, rho_eff_uc, 'Color', hex2rgb('#DEAA79'),... % [0.4940 0.1840 0.5560]
    'LineStyle', '-.', 'LineWidth', 3);
lgdtxt1{end} = 'Olfert \& Rogak (2019)';
hold on

% initialize figure for classified distribution shape vs. AAC setpoint 
f3 = figure(3);
f3.Position = [100, 100, 700, 800];
set(f3, 'color', 'white');
t3 = tiledlayout(2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

for i = ii
    
    % import saved workspace data for each condition
    fadd_wsp = strcat(fdir{i}, '\', fname{i}, '.mat');
    load(fadd_wsp, varnames{:})
    
    % double check distribution metrics
    dist_odias = adjust_dist(dist_odias);

    for j = 1 : size(dist_odias,1)

        % calculate effective density for each setpoint
        dist_odias(j).rho_eff = DAT.RHO_EFF(dist_odias(j).d_mode,...
            dat(j).da) ./ coef_rho(dist_odias(j).d_mode);
        
        % store aerodynamic setpoint
        dist_odias(j).da = dat(j).da;
    
    end

    % store the processed data
    dists{i} = dist_odias;
        
    % scatter plots of effective density vs. mobility diameter at each day
    figure(f1)
    nexttile(1)
    plt1{i} = scatter(cat(1,dists{i}.d_mode), cat(1,dists{i}.rho_eff),...
        mrksz(i), hex2rgb(clr{i}), mrk{i}, 'LineWidth', 1.5);
    hold on
    
    % plot scatters of GSD of mobility size distribution vs. mobility...
        % ...diameter at each day
    nexttile(2) % geometric sta
    scatter(cat(1,dists{i}.d_mode), cat(1,dists{i}.sigma_g),...
        mrksz(i), hex2rgb(clr{i}), mrk{i}, 'LineWidth', 1.5);
    hold on
    
    % make a description for the dataset
    lgdtxt1{i} = strcat(test_date(i), ',', {' '}, test_condition0(i));

    % plot the second to fourth moments of tandem distributions
    figure(f3)

    nexttile(1) % mode
    plt3{i} = scatter(cat(1,dists{i}.da), cat(1,dists{i}.d_mode),...
        mrksz(i), hex2rgb(clr{i}), mrk{i}, 'LineWidth', 1.5);
    hold on
    if i == ii(end)
        box on
        set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12,...
            'TickLength', [0.02 0.02], 'XScale', 'log', 'YScale', 'log')
        xlabel('$d_\mathrm{a}$ [nm]', 'interpreter', 'latex',...
            'FontSize', 16)
        ylabel('$d_\mathrm{mod} [-]$',...
            'interpreter', 'latex', 'FontSize', 16)
    end
    
    nexttile(2) % geometric standard deviation
    plt3{i} = scatter(cat(1,dists{i}.da), cat(1,dists{i}.sigma_g),...
        mrksz(i), hex2rgb(clr{i}), mrk{i}, 'LineWidth', 1.5);
    hold on
    if i == ii(end)
        box on
        set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12,...
            'TickLength', [0.02 0.02], 'XScale', 'log')
        xlabel('$d_\mathrm{a}$ [nm]', 'interpreter', 'latex',...
            'FontSize', 16)
        ylabel('$\sigma_\mathrm{m} [-]$',...
            'interpreter', 'latex', 'FontSize', 16)
    end    
    
    nexttile(3) % skewness
    scatter(cat(1,dists{i}.da), cat(1,dists{i}.skw),...
        mrksz(i), hex2rgb(clr{i}), mrk{i}, 'LineWidth', 1.5);
    hold on
    if i == ii(end)
        box on
        set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12,...
            'TickLength', [0.02 0.02], 'XScale', 'log')
        xlabel('$d_\mathrm{a}$ [nm]', 'interpreter', 'latex',...
            'FontSize', 16)
        ylabel('Skewness [-]', 'interpreter', 'latex', 'FontSize', 16)
    end    
    
    nexttile(4) % kurtosis
    scatter(cat(1,dists{i}.da), cat(1,dists{i}.krts),...
        mrksz(i), hex2rgb(clr{i}), mrk{i}, 'LineWidth', 1.5);
    hold on
    if i == ii(end)
        box on
        set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12,...
            'TickLength', [0.02 0.02], 'XScale', 'log')
        xlabel('$d_\mathrm{a}$ [nm]', 'interpreter', 'latex',...
            'FontSize', 16)
        ylabel('Kurtosis [-]', 'interpreter', 'latex', 'FontSize', 16)
    end    
    
    % clear redundant variables
    clear dist_odias dat rho_eff

end

% set apprearances for figure 1
figure(f1)
nexttile(1)
box on
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12,...
    'TickLength', [0.02 0.02], 'XScale', 'log', 'YScale', 'log')
xlabel('$d_\mathrm{m}$ [nm]', 'interpreter', 'latex',...
    'FontSize', 16)
ylabel('$\rho_\mathrm{eff} [\mathrm{kg}/\mathrm{m}^3]$',...
    'interpreter', 'latex', 'FontSize', 16)
hold on
nexttile(2)
box on
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12,...
    'TickLength', [0.02 0.02], 'XScale', 'log')
xlabel('$d_\mathrm{m}$ [nm]', 'interpreter', 'latex',...
    'FontSize', 16)
ylabel('$\sigma_\mathrm{m} [-]$', 'interpreter', 'latex', 'FontSize', 16)

%% present data in groups if requested by user

[~, idx] = unique(grp, 'stable');
test_condition = test_condition0(idx);

if exist('ind0_grp', 'var') && ~isempty(ind0_grp)

    % initialize the figure for da vs. dm curve fitting
    f4 = figure(4);
    f4.Position = [400, 100, 500, 500];
    set(f4, 'color', 'white');
    plt4 = cell(n_grp,1);
    
    % placeholder for curve fitting data
    fit4 = cell(n_grp,1);
    da_fit = cell(n_grp,1);
    dm_fit = cell(n_grp,1);

    iii = []; % initialize index for legends of figure 1

    % initialize the figure for ensemble groups
    f2 = figure(2);
    f2.Position = [150, 150, 900, 500];
    set(f2, 'color', 'white');
    t2 = tiledlayout(1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
    
    % initialize plot and legend placeholders
    plt21 = cell(n_grp + 2, 1);
    lgdtxt2 = cell(n_grp + 2, 1);
    
    nexttile(1)
    % plot universal correlation
    plt21{end-1} = plot(dm_uc, rho_eff_uc, 'Color', hex2rgb('#DEAA79'),... % [0.4940 0.1840 0.5560]
        'LineStyle', '-.', 'LineWidth', 3);
    lgdtxt2{end-1} = 'Olfert \& Rogak (2019)';
    hold on

    %%% Sipkens and Corbin (2024) correlation for collapsed soot
    % correlation constants
    rho_eff_100_ca = 786;    % kg/m^3
    rho_eff_c = 651;         % kg/m^3
    zeta_ca = 2.44;
    dm_transition = 140;     % nm
    dm_100 = 100;            % nm (reference point)

    % Preallocate output
    rho_eff_col = zeros(size(dm_uc));

    % Apply piecewise definition
    rho_eff_col(dm_uc < dm_transition) = rho_eff_100_ca .*...
        (dm_uc(dm_uc < dm_transition) ./ dm_100) .^ (zeta_ca - 3);
    rho_eff_col(dm_uc >= dm_transition) = rho_eff_c;

    % plot correlation
    plt21{end} = plot(dm_uc, rho_eff_col, 'Color', hex2rgb('#9FB3DF'),...
        'LineStyle', ':', 'LineWidth', 2.5);
    lgdtxt2{end} = 'Sipkens \& Corbin (2024)';
    %%%
    
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
            dist_grp(k).sigma_g = [dist_grp(k).sigma_g,...
                dists{ind_grp{k}(kk)}.sigma_g];
            dist_grp(k).da = [dist_grp(k).da,...
                dists{ind_grp{k}(kk)}.da];

            % sort elements based on the order of da
            [dist_grp(k).da, id_sort] = sort(dist_grp(k).da);
            dist_grp(k).d_mode = dist_grp(k).d_mode(id_sort);
            dist_grp(k).d_gm = dist_grp(k).d_gm(id_sort);
            dist_grp(k).rho_eff = dist_grp(k).rho_eff(id_sort);
            dist_grp(k).sigma_g = dist_grp(k).sigma_g(id_sort);
            
            % sort out indices for order of figure 1 legend
            iii = [iii, ind_grp{k}(kk)];
            
        end
                      
        % plot the grouped effective density datasets
        figure(f2)
        nexttile(1)
        plt21{k} = scatter(dist_grp(k).d_mode, dist_grp(k).rho_eff,...
            mrksz2(k), hex2rgb(clr2{k}{1}), mrk2{k}, 'LineWidth', 1.5);
        hold on

        % lgdtxt2{k} = strcat('Group', {' '}, num2str(k));
        lgdtxt2{k} = test_condition(k);

        % set apprearances
        box on
        set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 16,...
            'TickLength', [0.02 0.02], 'XScale', 'log', 'YScale', 'log')
        xlabel('$d_\mathrm{m}$ [nm]', 'interpreter', 'latex',...
            'FontSize', 24)
        ylabel('$\rho_\mathrm{eff} [\mathrm{kg}/\mathrm{m}^3]$',...
            'interpreter', 'latex', 'FontSize', 24)
        xlim([0.8 * min(cat(2,dist_grp.d_mode))...
            1.2 * max(cat(2,dist_grp.d_mode))])
        ylim([0.8 * min(cat(2,dist_grp.rho_eff))...
            1.5 * max(cat(2,dist_grp.rho_eff))])

        nexttile(2)
        plt22{k} = scatter(dist_grp(k).d_mode, dist_grp(k).sigma_g,...
            mrksz2(k), hex2rgb(clr2{k}{1}), mrk2{k}, 'LineWidth', 1.5);
        hold on

        % set apprearances
        box on
        set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 16,...
            'TickLength', [0.02 0.02], 'XScale', 'log')
        xlabel('$d_\mathrm{m}$ [nm]', 'interpreter', 'latex',...
            'FontSize', 24)
        ylabel('$\sigma_\mathrm{m}$  [-]', 'interpreter', 'latex',...
            'FontSize', 24)
        
        figure(f4) % curve fit to dm vs da data

        % % use spline
        % t = linspace(0,1,length(cat(1,dist_grp(k).da)));
        % p = 0.001; % Smoothing parameter (adjust as needed)
        % sp = spaps(t, [dist_grp(k).da; dist_grp(k).d_mode], p);
        % tq = linspace(0,1,100);
        % fit4{k} = fnval(sp, tq);
        % plt4 {k} = plot(fit4{k}(1,:), fit4{k}(2,:), 'Color',...
        %     hex2rgb(clr2{k}{2}), 'LineStyle', '-', 'LineWidth', 1.5);
        % hold on

        % use optimized polyfit model
        fit4{k} = fit(dist_grp(k).da', dist_grp(k).d_gm', fittype{k},...
            'Normalize', 'on', 'Robust', 'on');
        da_fit{k} = linspace(min(dist_grp(k).da),...
            max(dist_grp(k).da), 100)'; % set extrapolated range
        dm_fit{k} = fit4{k}(da_fit{k}); % evaluate model on new points
        plt4 {k} = plot(da_fit{k}, dm_fit{k}, 'Color',...
            hex2rgb(clr2{k}{2}), 'LineStyle', '-', 'LineWidth', 1.5);
        hold on

        % plot original data
        scatter(dist_grp(k).da, dist_grp(k).d_mode,...
            mrksz(k), hex2rgb(clr2{k}{1}), mrk{k}, 'LineWidth', 1.5)
        
    end
    
    lgd2 = legend(cat(2, plt21{:}), cat(2, lgdtxt2{:}), 'interpreter',...
        'latex', 'FontSize', 20, 'Location', 'northoutside',...
        'NumColumns', 3);
    lgd2.Layout.Tile = 'south';
    
    % adjust the bounds in non-grouped effective density figure
    figure(f1)
    nexttile(1)
    xlim([0.8 * min(cat(2,dist_grp.d_mode))...
        1.2 * max(cat(2,dist_grp.d_mode))])
    ylim([0.8 * min(cat(2,dist_grp.rho_eff))...
        1.2 * max(cat(2,dist_grp.rho_eff))])

    % sort out order of legends in non-grouped figures
    lgd1 = legend(cat(2, plt1{iii}, plt1{end}), cat(2, lgdtxt1{iii},...
    lgdtxt1{end}), 'interpreter', 'latex', 'FontSize', 12);
    lgd1.Layout.Tile = 'south';
    lgd1.NumColumns = 2;
    
    figure(f3)
    lgd3 = legend(cat(2, plt3{iii}), cat(2, lgdtxt1{iii}),...
    'interpreter', 'latex', 'FontSize', 12);
    lgd3.Layout.Tile = 'south';
    lgd3.NumColumns = 2;

    figure(f4)
    box on
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12,...
        'TickLength', [0.02 0.02], 'XScale', 'log', 'YScale', 'log')
    xlabel('$d_\mathrm{a}$ [nm]', 'interpreter', 'latex',...
        'FontSize', 16)
    ylabel('$d_\mathrm{m} [nm]$', 'interpreter', 'latex', 'FontSize', 16)
    lgd2 = legend(cat(2, plt4{:}), cat(2, lgdtxt2{1:end-1}),...
        'interpreter', 'latex', 'FontSize', 12, 'Location', 'southeast');    
    
else

    % just print legend for the non-grouped figures
    legend(cat(2, plt1{:}), cat(2, lgdtxt1{:}), 'interpreter', 'latex',...
        'FontSize', 12, 'Location', 'southoutside', 'NumColumns', 3);

end

%% save data

% make directory for results to be saved
dir_out_name = 'Effective-Density-Compiled';
fdir_out = strcat('outputs\', dir_out_name);
if ~isfolder(fdir_out); mkdir(fdir_out); end

% customize the name of directory
fname_out = strcat(dir_out_name, '_', datestr(datetime('now'))); %#ok<DATST>
fname_out = regexprep(fname_out, ':', '-');
fname_out = regexprep(fname_out, ' ', '_');

% save MATLAB worspace
% save(strcat(fdir_out, '\', fname_out, '.mat'));

function dist_odias = adjust_dist(dist_odias)
% recalculate parameters of size distribution

for i = 1 : length(dist_odias)

    % remove negative counts (artifacts)
    iii = dist_odias(i).x < 0;
    dist_odias(i).x2 = dist_odias(i).x;
    dist_odias(i).d2 = dist_odias(i).d;
    if nnz(iii)    
        dist_odias{i}.x2(iii) = [];
        dist_odias{i}.d2(iii) = [];
    end

    % find local modes of distribution and remove noise
    [dn_dlogd_mode, dist_odias(i).d_mode] = findpeaks(dist_odias(i).x2,...
        dist_odias(i).d2); 
    dist_odias(i).d_mode(dn_dlogd_mode / max(dn_dlogd_mode) < 0.1) = [];

    % find geometric mean (GM) and geometric standard deviation (GSD)
    w = dist_odias(i).x2 / sum(dist_odias(i).x2); % normalize dn/dlog(d) to get weights
    dist_odias(i).d_gm = 10^(sum(w .* log10(dist_odias(i).d2))); % GM
    dist_odias(i).sigma_g = 10^(sqrt(sum(w .* (log10(dist_odias(i).d2) -...
        log10(dist_odias(i).d_gm)).^2))); % GSD

    % total concentration (i.e. area below the size distribution curve)
    dist_odias(i).n_tot = trapz(log10(dist_odias(i).d2), dist_odias(i).x2);
    
    % skewness in log-space (0 for a normal distribution, positive for...
    % ...right-skewed, negative for left-skewed)
    dist_odias(i).skw = sum(w .* (log10(dist_odias(i).d2) -...
        log10(dist_odias(i).d_gm)).^3) / ((log10(dist_odias(i).sigma_g))^3);
    
    % kurtosis in log-space (3 for a normal distribution, > 3 for...
    % ...heavy tails, < 3 for light tails)
    dist_odias(i).krts = sum(w .* (log10(dist_odias(i).d2) -...
        log10(dist_odias(i).d_gm)).^4) / ((log10(dist_odias(i).sigma_g))^4);
    
    % find the highest peak among the modes at each AAC setpoint
    dist_odias = hn.selectpeak(dist_odias);
    
end

end
