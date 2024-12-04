
clear;
close all;
clc;
warning('off')

addpath cmap;

%%% read input parameters from the user %%%
tname0 = '13SEP24'; % name of excel datasheet for cases to be imported
f_in_add = 'inputs\odias_params_new.xlsx'; % address of folder containing main user info
opts = detectImportOptions(f_in_add);
opts = setvartype(opts, 'char'); % set excel file import option
intab = readtable(f_in_add, 'UseExcel', true, 'Sheet', tname0); % read the excel file
fname = intab.fname{1}; % identify name of smps files that contain raw particle counts
tname = intab.tname{1}; % name of datasheet that the above smps data are located in 
fname_lbl = intab.fname_lbl{1}; % label for the data to be imported
n_dat = intab.n_dat(1); % number of test cases
lambda_tk1 = intab.lambda_tk1; % a parameter that controls smoothing level in Tikhonov method
clr = intab.clr; % assign color for plot lines
linstl = intab.stl; % assign line style for plots
d_lim_1 = intab.d_lim_1; % lower limit on mobility diameter (for inversion)
d_lim_2 = intab.d_lim_2; % upper limit on mobility diameter
dm_rowid = intab.dm_rowid; % row id in SMPS file that corresponds to mobility setpoints
%%% %%%

dat = DAT(); % initialize the class of SMPS data

% generate a text description for the test
rigstr = DAT.LABELMAKER(fname_lbl, tname);

prop = tfer.prop_dma; % load default SMPS properties

% make placeholders to store mobility setpoints and counts for use in log-scale
dm2 = cell(n_dat, 1);
dN2 = cell(n_dat, 1);

% initialize mobility size distribution figure for all cases
f1 = figure(1);
f1.Position = [100, 100, 1000, 600];
set(f1, 'color', 'white');
t1 = tiledlayout(1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

plt1 = cell(n_dat,1); % placeholder to store plots for assigning legends

tools.textbar([0, n_dat]); % initialize progress textbar

for i = 1 : n_dat

    nexttile(1) % non-log scales
    
    % import SMPS files and average particle counts
    dat(i) = DAT(i, fname, tname, dm_rowid(i));

    % adjust the description text
    % rigstr{i} = DAT.LABELADJ(rigstr{i});
    
    % reconstruct mobility points
    d = logspace(log10(d_lim_1(i)), log10(d_lim_2(i)), 500)';
    
    % correct sheath flows
    prop.Q_c = dat(i).Qsh_dma;
    prop.Q_m = dat(i).Qsh_dma;

    % build inversion kernel
    A = kernel.gen_smps((dat(i).dm)', d, [], prop);
    
    % assign mean and std of particle counts
    b = (dat(i).dN)';
    Lb = bsxfun(@rdivide, eye(length(b)), (dat(i).sigma)');

    disp(' ');

    % perform inversion using 1st order Tikhonov
    % disp('Running Tikhonov (1st) ...');
    [x_tk1, ~, ~, Gpo_inv_tk1] = ...
        invert.tikhonov(Lb * A, Lb * b, lambda_tk1(i), 1);
    Gpo_tk1 = inv(Gpo_inv_tk1);
    % tools.textdone();
    % disp(' ');

    % display the inverted mobility size distribution based on Tikhonov
    plt1{i} =  tools.plotci_HN(d, x_tk1, Gpo_tk1, [], hex2rgb(clr{i}), linstl{i});
    hold on
    
    % rigstr{i} = strrep(rigstr{i},'nof','n/f'); % minor label adjustment

    % apply the axis configurations
    if i == n_dat
        set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12,...
            'TickLength', [0.02 0.02], 'XScale', 'log')    
        xlim([15,1000])
        xlabel('$d_\mathrm{m}$ [nm]', 'interpreter', 'latex',...
            'FontSize', 16)
        ylabel('$\mathrm{d}n/\mathrm{dlog}(d_\mathrm{m})$',...
            'interpreter', 'latex', 'FontSize', 16)
    end

    nexttile(2) % log scales

    % remove zero counts for plotting in log scale
    dm2{i} = d(x_tk1 > 0);
    dN2{i} = x_tk1(x_tk1 > 0);
    
    % plot nonzero scales for display in log-log scale
    opts2.log = 'on';
    tools.plotci_HN(dm2{i}, dN2{i}, Gpo_tk1, [], hex2rgb(clr{i}),...
        linstl{i}, opts2)
    hold on

    % apply the log-log scale and other axis configurations
    if i == n_dat
        set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12,...
            'TickLength', [0.02 0.02], 'XScale', 'log', 'YScale', 'log')    
        xlim([15,1000])
        ylim([1e2 1.2 * max(ylim)])
        xlabel('$d_\mathrm{m}$ [nm]', 'interpreter', 'latex',...
            'FontSize', 16)
    end

    tools.textbar([i, n_dat]); % update textbar
end

lgd1 = legend(cat(2, plt1{:}), cat(2, rigstr(:)), 'interpreter', 'latex',...
    'FontSize', 12);
lgd1.Layout.Tile = 'south';
lgd1.NumColumns = 2;


