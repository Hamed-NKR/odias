
clear;
close all;
clc;
warning('off')

addpath cmap;

%%% read user input parameters %%%
tname0 = '26&30JUL24'; % datasheet to be imported
f_in_add = 'inputs\odias_params.xlsx'; 
opts = detectImportOptions(f_in_add);
opts = setvartype(opts, 'char');
intab = readtable(f_in_add, 'UseExcel', true, 'Sheet', tname0);
fname = intab.fname{1}; % filename to be imported
tname = intab.tname{1}; % tabname that contains data of interest
fname_lbl = intab.fname_lbl{1};
n_dat = intab.n_dat(1);
lambda_tk1 = intab.lambda_tk1;
clr = intab.clr;
d_lim_1 = intab.d_lim_1;
d_lim_2 = intab.d_lim_2;
%%% %%%


dat = DAT(); % initialize the class of DMA data

% generate a text description for the test
rigstr = DAT.LABELMAKER(fname_lbl, tname);

prop = tfer.prop_dma; % load default dma properties

% initialize figure
figure;
h = gcf;
h.Position = [0, 0, ceil(n_dat/8) * 650, 450];
set(h, 'color', 'white');

plt = cell(n_dat,1); % placeholder to store plots for assigning legends

tools.textbar([0, n_dat]); % initialize textbar

for i = 1 : n_dat
    
    % import DMA files and average particle counts
    dat(i) = DAT(i, fname, tname);

    % adjust the description text
    rigstr{i} = DAT.LABELADJ(rigstr{i});
    
    % reconstruction points
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

    % inversion using 1st order Tikhonov
    % disp('Running Tikhonov (1st) ...');
    [x_tk1, ~, ~, Gpo_inv_tk1] = ...
        invert.tikhonov(Lb * A, Lb * b, lambda_tk1(i), 1);
    Gpo_tk1 = inv(Gpo_inv_tk1);
    % tools.textdone();
    % disp(' ');
    
    % plot mobility size distribution
    plt{i} =  tools.plotci(d, x_tk1, Gpo_tk1, [], hex2rgb(clr{i}));
    hold on
    
    rigstr{i} = strrep(rigstr{i},'nof','n/f'); % minor label adjustment

    tools.textbar([i, n_dat]); % update textbar
end

% plot appearance configs
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12,...
    'TickLength', [0.02 0.02])
xlim([15,1500])
xlabel('$d_\mathrm{m}$ [nm]', 'interpreter', 'latex',...
    'FontSize', 16)
ylabel('$\mathrm{d}n/\mathrm{dlog}(d_\mathrm{m})$',...
    'interpreter', 'latex', 'FontSize', 16)
% legend(cat(2, plt{:}), cat(2, rigstr{:}), 'interpreter', 'latex',...
%     'Location', 'eastoutside', 'NumColumns', ceil(n_dat/8), 'FontSize', 12)


