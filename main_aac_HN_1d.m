clear;
close all;
clc;
warning('off')

%% user inputs

% name and address of batch file
fadd_batch = 'inputs\aac_params_new.xlsx';

%% read batch file

% set data type and import excel file
opts_read = detectImportOptions(fadd_batch);
opts_read = setvartype(opts_read, 'char');
intab = readtable(fadd_batch, 'UseExcel', true);

% sort out batch file parameters
n0_dat = intab.n_dat(1); % total number of data given by the user
ii0 = 1 : n0_dat; % indices of the input data
fdir = intab.fdir(ii0); % folder addresses of aac files
fname = intab.fname(ii0); % filenames of aac files
clr1 = intab.clr1(ii0); % assign colors to lines
clr2 = intab.clr2(ii0); % assign colors to plot shades
linstl = intab.linstl(ii0); % assign styles for plot lines
linsz = intab.linsz(ii0); % assign line thicknesses
test_date = intab.date(ii0); % date of measurement
test_condition = intab.condt(ii0); % condition of measurement
DR = intab.DR(ii0); % dilution ratio
incld = intab.incld(ii0); % determine whether a dataset would be...
    % ...post-processed
ii = ii0(incld); % indices of data selected by the user for post-processing
n_dat = length(ii); % number of datasets to be analyzed

%% Load and plot AAC distribution data

% initialize figure
f1 = figure(1);
f1.Position = [50, 50, 500, 500];
set(f1, 'color', 'white');

tab_aac = cell(n_dat,1); % placeholder for raw aac files

% placeholders for locations of size distribution data in aac files
start_idx = cell(n_dat,1);
end_idx = cell(n_dat,1);

% placeholders for data from individual distributions
dists_aac = cell(n_dat,1);

% placeholder for distribution plots
plt = cell(n_dat,1);

for i = ii
    
    % make the full file address
    fadd_aac = strcat(fdir{i}, '\', fname{i}, '.csv');
    
    % set configs for import
    opts = detectImportOptions(fadd_aac);
    opts = setvartype(opts, 'double');

    % import aac files
    tab_aac{i} = readtable(fadd_aac, opts);

    % identify the start and end points for data of each distribution
    is_non_NaN = ~isnan(table2array(tab_aac{i}(:,3))); % identify non-NaN...
        % ...indices
    start_idx{i} = find(diff([0; is_non_NaN]) == 1); % find the start indices...
        % ...of non-NaN series
    end_idx{i} = find(diff([is_non_NaN; 0]) == -1); % find the end indices...
        % ...of non-NaN series
    
    % number of scans existing in the aac file
    dists_aac{i}.n_scan = length(start_idx{i});

    % initialize placeholder for aerodynamic setpoints and corresponding...
        % ...particle counts
    dists_aac{i}.da = cell(dists_aac{i}.n_scan,1);
    dists_aac{i}.dn_dlogda = cell(dists_aac{i}.n_scan,1);  
    
    for j = 1 : dists_aac{i}.n_scan

        dists_aac{i}.da0{j} = table2array(tab_aac{i}(start_idx{i}(j):...
            end_idx{i}(j),3));
        dists_aac{i}.dn0_dlogda{j} = DR(i) * table2array(tab_aac{i}(start_idx{i}(j):...
            end_idx{i}(j),5));
        
    end
    
    [dists_aac{i}.dn_dlogda, dists_aac{i}.da, dists_aac{i}.ci_dn_dlogd] = ...
        hn.avg_dist_2(dists_aac{i}.da0, dists_aac{i}.dn0_dlogda);
    ub_dn_dlogda = dists_aac{i}.dn_dlogda + dists_aac{i}.ci_dn_dlogd;
    lb_dn_dlogda = dists_aac{i}.dn_dlogda - dists_aac{i}.ci_dn_dlogd;
    
    iii = dists_aac{i}.dn_dlogda < 0;
    if nnz(iii)
        dists_aac{i}.da(iii) = [];
        dists_aac{i}.dn_dlogda(iii) = [];
        dists_aac{i}.ci_dn_dlogda(iii) = [];
        dists_aac{i}.ub_dn_dlogda(iii) = [];
        dists_aac{i}.lb_dn_dlogda(iii) = [];
    end

    plt{i} = plot(dists_aac{i}.da, dists_aac{i}.dn_dlogda, 'Color',...
        hex2rgb(clr1{i}), 'LineStyle', linstl{i}, 'LineWidth', linsz(i));
    hold on
    fill([dists_aac{i}.da, fliplr(dists_aac{i}.da)],...
        [ub_dn_dlogda, fliplr(lb_dn_dlogda)],...
        hex2rgb(clr2{i}), 'EdgeColor', 'none', 'FaceAlpha', 0.5);

    % total concentration
    dists_aac{i}.n_tot = trapz(log10(dists_aac{i}.da),...
        dists_aac{i}.dn_dlogda);
    
    % weights for calculating geometric mean
    w = dists_aac{i}.dn_dlogda / sum(dists_aac{i}.dn_dlogda);
    
    % geometric mean
    dists_aac{i}.gm_da = 10^(sum(w .* log10(dists_aac{i}.da)));
    
    % geometric standard deviation
    dists_aac{i}.gsd_da = 10^(sqrt(sum(w .* (log10(dists_aac{i}.da) -...
        log10(dists_aac{i}.gm_da)).^2)));    
    
end

set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12,...
    'TickLength', [0.015 0.015], 'XScale', 'log')
xlim([35,600])
ylim([-Inf, 1.1e8])
xlabel('$d_\mathrm{ae}$ [nm]', 'interpreter', 'latex',...
    'FontSize', 16)
ylabel('$\mathrm{d}n/\mathrm{dlog}(d_\mathrm{ae}) [\#/\mathrm{cm}^3]$',...
    'interpreter', 'latex', 'FontSize', 16)
legend(cat(2, plt{:}), cat(2, test_condition(:)), 'interpreter', 'latex',...
    'FontSize', 12, 'Location', 'northeast');

